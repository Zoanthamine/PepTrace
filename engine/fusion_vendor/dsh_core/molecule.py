"""分子/中间体模型（T04）。

执行计划第 3 节 T04 要求建立：链、端基、修饰连接、保护态、**固定电荷、同位素、起始·目标·中间体**模型，
并扩充组件约束。state.py 已负责「链 + 端基 + 修饰 + 保护态 → 规范态 → 分子式」；
本模块补上另外三件事：

1. **起始 / 中间体 / 终产物**统一成 `Molecule` 节点：起始肽、加一个事件后的中间体、再叠一个事件的终产物，
   都是同一种对象，只有这样组合搜索（T06）才能按「当前分子状态」继续反应，而不是把净质量差简单相加。
2. **组件四类分开用**：游离氨基酸、肽内残基、试剂本体、共价净增量。组件库里四张表都有，
   但此前计算路径只用了「残基 + H2O」；游离氨基酸（如水解出的单个氨基酸）与试剂本体（如哌啶、
   TFA 本身进入产物的情形）需要独立取用，不能拿残基质量当游离氨基酸质量。
3. **质量口径齐全**：单同位素质量、平均质量、同位素包络、离子 m/z（含固定电荷），一处取齐。

重要边界
--------
`Molecule.mass_only=True` 表示只知道质量、不知道元素组成（对应规则库里 156 条
`mass_delta_value` 型规则：只给数值不给净式）。这类节点能做质量累加与匹配，但**不能**算平均质量、
同位素包络，也不能做固定电荷换算，调用这些功能会抛 `MoleculeError`，不会静默给个错数。
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .components import ComponentLibrary
from .formula import Formula, delta_str
from .ions import FIXED_CHARGE_GROUPS, IonSpecies
from .isotopes import average_mass_of, envelope


class MoleculeError(ValueError):
    """分子节点操作不合法（如对无组成节点求平均质量）。"""


def parse_delta_string(text: str) -> Formula:
    """把 "C+1O+2" / "H-2" / "H2S-1" / "H2O" / "-H2O" 这类净式写成 Formula。

    支持两种写法：
    - 逐元素带符号：`C+1O+2`、`H-2Br-2`、`H2S-1`（空格忽略）；
    - 整体正配平：`H2O`；或整体取负：`-H2O`。
    认不出就抛错，不做猜测。
    """
    import re

    s = (text or "").strip().replace(" ", "")
    if not s:
        raise MoleculeError("净式为空")
    outer = 1
    if s[0] in "+-":
        outer = -1 if s[0] == "-" else 1
        s = s[1:]

    token = re.compile(r"([A-Z][a-z]?)(\d*)([+-]\d+)?")
    counts: dict[str, int] = {}
    pos = 0
    while pos < len(s):
        m = token.match(s, pos)
        if not m:
            raise MoleculeError(f"无法解析净式 {text!r}（卡在位置 {pos}：{s[pos:pos + 6]!r}）")
        el, num, signed = m.group(1), m.group(2), m.group(3)
        n = int(signed) if signed else outer * (int(num) if num else 1)
        counts[el] = counts.get(el, 0) + n
        pos = m.end()
    if not counts or all(v == 0 for v in counts.values()):
        raise MoleculeError(f"净式 {text!r} 解析结果为空")
    try:
        return Formula.of({k: v for k, v in counts.items() if v})
    except Exception as exc:  # 未知元素等统一收敛为 MoleculeError
        raise MoleculeError(f"净式 {text!r} 含非法元素或下标：{exc}") from exc


@dataclass(frozen=True)
class Molecule:
    """反应网络里的一个分子节点。

    stage: "start"（起始物）| "intermediate"（中间体）| "target"（观测到的目标质量）| "final"（终产物）
    """

    label: str
    mass: float                                  # 中性单同位素质量（Da）
    formula: Formula | None = None
    stage: str = "start"
    sequence: str | None = None
    fixed_charge: int = 0
    edits: tuple[str, ...] = ()                  # 到达该节点所做的变换说明
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------ 构造
    @classmethod
    def from_state(cls, state, lib: ComponentLibrary, label: str = "起始肽") -> "Molecule":
        f = state.formula(lib)
        return cls(label=label, mass=f.mass(), formula=f, stage="start",
                   sequence=state.sequence,
                   meta={"state_describe": state.describe()})

    @classmethod
    def from_free_amino_acid(cls, res: str, lib: ComponentLibrary) -> "Molecule":
        """游离氨基酸（不是肽内残基！两者差一个 H2O）。"""
        try:
            f = lib.free_amino_acids[res.upper()].formula
        except KeyError as exc:
            raise MoleculeError(f"组件库没有游离氨基酸 {res!r}") from exc
        return cls(label=f"游离氨基酸 {res.upper()}", mass=f.mass(), formula=f,
                   stage="final", meta={"component_kind": "free_amino_acid"})

    @classmethod
    def from_reagent(cls, name: str, lib: ComponentLibrary) -> "Molecule":
        """试剂本体（中性分子），用于"试剂整体加成/进入产物"一类事件。"""
        key = name
        if key not in lib.neutral:
            low = {k.lower(): k for k in lib.neutral}
            if name.lower() in low:
                key = low[name.lower()]
            else:
                raise MoleculeError(f"组件库没有中性试剂 {name!r}；可选 {sorted(lib.neutral)[:8]}…")
        f = lib.neutral[key].formula
        return cls(label=f"试剂 {key}", mass=f.mass(), formula=f, stage="start",
                   meta={"component_kind": "reagent"})

    @classmethod
    def from_sequence_all_residues(cls, sequence: str, lib: ComponentLibrary) -> "Molecule":
        """把序列当"游离氨基酸之和"（即完全水解产物），用于判断水解类观测质量。"""
        total = Formula.of({})
        for ch in sequence.upper():
            if ch not in lib.free_amino_acids:
                raise MoleculeError(f"组件库没有游离氨基酸 {ch!r}")
            total = total + lib.free_amino_acids[ch].formula
        return cls(label=f"完全水解 {sequence.upper()}", mass=total.mass(), formula=total,
                   stage="final", sequence=sequence.upper(),
                   meta={"component_kind": "free_amino_acids_sum"})

    # ------------------------------------------------------------ 变换
    def apply_delta(self, delta_formula: Formula | dict | str | None,
                    delta_mass: float | None = None,
                    label: str | None = None,
                    edit: str | None = None,
                    stage: str = "intermediate") -> "Molecule":
        """按净式（优先）或净质量差生成下一节点。

        净式优先：只有净式能保证元素组成可继续传递；只给数值时节点降级为 mass_only。
        """
        delta: Formula | None = None
        if delta_formula is not None:
            if isinstance(delta_formula, Formula):
                delta = delta_formula
            elif isinstance(delta_formula, dict):
                delta = Formula.of(delta_formula)
            else:
                delta = parse_delta_string(str(delta_formula))

        if delta is not None:
            if self.formula is None:
                # 已知净式但当前节点无组成：质量可加，组成仍未知
                return Molecule(label=label or f"{self.label} {delta_str(delta.as_dict())}",
                                mass=self.mass + delta.mass(), formula=None,
                                stage=stage, sequence=self.sequence,
                                fixed_charge=self.fixed_charge,
                                edits=self.edits + ((edit,) if edit else ()),
                                meta={**self.meta, "mass_only": True})
            new_formula = self.formula + delta
            return Molecule(label=label or f"{self.label} {delta_str(delta.as_dict())}",
                            mass=new_formula.mass(), formula=new_formula, stage=stage,
                            sequence=self.sequence, fixed_charge=self.fixed_charge,
                            edits=self.edits + ((edit,) if edit else ()),
                            meta=dict(self.meta))

        if delta_mass is None:
            raise MoleculeError("apply_delta 需要净式或净质量差之一")
        return Molecule(label=label or f"{self.label} {delta_mass:+.6f} Da",
                        mass=self.mass + float(delta_mass), formula=None, stage=stage,
                        sequence=self.sequence, fixed_charge=self.fixed_charge,
                        edits=self.edits + ((edit,) if edit else ()),
                        meta={**self.meta, "mass_only": True})

    def combine(self, other: "Molecule", label: str, edit: str | None = None) -> "Molecule":
        """两个节点合并（如肽 + 试剂整体加成）。只有双方都有组成时结果才有组成。"""
        formula = (self.formula + other.formula
                   if self.formula is not None and other.formula is not None else None)
        return Molecule(label=label, mass=self.mass + other.mass, formula=formula,
                        stage="intermediate", sequence=self.sequence,
                        fixed_charge=self.fixed_charge + other.fixed_charge,
                        edits=self.edits + ((edit,) if edit else ()),
                        meta={**self.meta, "mass_only": formula is None})

    # ------------------------------------------------------------ 质量口径
    @property
    def mass_only(self) -> bool:
        return self.formula is None

    def average_mass(self) -> float:
        if self.formula is None:
            raise MoleculeError(f"{self.label} 只有质量没有组成，无法计算平均质量；"
                                "对应规则只给了数值型净质量差")
        return average_mass_of(self.formula)

    def isotope_envelope(self, max_offset: int = 2) -> list[dict]:
        if self.formula is None:
            raise MoleculeError(f"{self.label} 只有质量没有组成，无法计算同位素包络")
        return envelope(self.formula, max_offset)

    def ion(self, charge: int, adduct: str = "H",
            fixed_charge: int | None = None) -> IonSpecies:
        fq = self.fixed_charge if fixed_charge is None else fixed_charge
        if fq and self.formula is None:
            raise MoleculeError(f"{self.label} 只有质量没有组成，不能把固定电荷换算到 m/z")
        return IonSpecies(self.mass, charge, adduct, fq)

    def mz(self, charge: int, adduct: str = "H") -> float:
        return self.ion(charge, adduct).mz()

    def describe(self) -> str:
        bits = [self.label, f"{self.mass:.6f} Da"]
        if self.formula is not None:
            bits.append(self.formula.name())
        else:
            bits.append("仅质量/无组成")
        if self.fixed_charge:
            bits.append(f"固定电荷 {self.fixed_charge:+d}")
        if self.edits:
            bits.append("路径：" + " → ".join(self.edits))
        return "；".join(bits)


def fixed_charge_catalog() -> dict[str, dict]:
    """固定电荷基团目录（供输入/报告引用，未登记的一律不假设）。"""
    return dict(FIXED_CHARGE_GROUPS)


def apply_rule_delta(mol: Molecule, rule, stage: str = "intermediate") -> Molecule:
    """把一条编译后规则（rules/CHxx.json 的 dict 或 Rule 对象）作用到分子节点。"""
    get = (lambda k: rule.get(k)) if isinstance(rule, dict) else (lambda k: getattr(rule, k, None))
    delta_formula = get("delta_formula")
    delta_mass = get("delta_mass_mono_da")
    rid = get("rule_id")
    name = get("name_cn") or rid
    return mol.apply_delta(delta_formula, delta_mass, label=f"{mol.label} →{name}",
                           edit=f"{rid} {name}", stage=stage)