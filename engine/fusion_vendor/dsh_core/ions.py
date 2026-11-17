"""离子化与固定电荷模型（T04）。

执行计划第 6 节明确要求「固定电荷独立建模」。为什么不能只写 [M+nH]n+：
  * 季铵、锍、吡啶𬭩等**永久电荷**基团自带 z_fixed 个正电荷，离子总电荷 z_total = z_fixed + z_proton，
    但质量只加 z_proton 个质子——把两者混在一起算会让 m/z 系统偏大。
  * 负离子模式（[M−nH]n−）与加合离子（Na+、K+、NH4+）走的是不同加合项，不能共用一套公式。
  * 两个实测峰的差值 ∆m/z 必须按电荷换算回中性质量差，并传播误差；这一层单独提供函数。

质量约定
--------
- 所有质量都是中性单同位素质量（Da）。
- 质子质量用 PROTON（氢原子质量减电子），电子质量 ELECTRON。
- 误差方向统一「实测 − 理论」（执行计划第 6 节）。
"""
from __future__ import annotations

from dataclasses import dataclass

from .components import ComponentLibrary
from .formula import ELECTRON, PROTON

# 常见加合离子对应的中性质量增量（离子的质量，不含电荷补偿）
ADDUCTS: dict[str, dict] = {
    "H":    {"label": "[M+nH]ⁿ⁺",     "mass": PROTON,   "charge": +1, "protons": True},
    "Na":   {"label": "[M+nNa]ⁿ⁺",    "mass": 22.98976928 - ELECTRON, "charge": +1, "protons": False},
    "K":    {"label": "[M+nK]ⁿ⁺",     "mass": 38.96370649 - ELECTRON, "charge": +1, "protons": False},
    "NH4":  {"label": "[M+nNH4]ⁿ⁺",   "mass": 18.03382555 - ELECTRON, "charge": +1, "protons": False},
    "H-neg": {"label": "[M−nH]ⁿ⁻",    "mass": -PROTON,  "charge": -1, "protons": True},
}

# 离子化残基/端基的近似 pKa（只用于"可能带几个电荷"的提示，不用于定量 pKa 计算）
_PKA = {
    "Nterm": 9.6, "Cterm": 2.4, "K": 10.5, "R": 12.5, "H": 6.0,
    "D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1,
}

# 固定电荷基团：名称 → (固定电荷数, 说明)。只放化学上确实带永久电荷的基团。
FIXED_CHARGE_GROUPS: dict[str, dict] = {
    "QuaternaryAmmonium": {"charge": +1, "note": "季铵（如 N,N,N-三甲基铵）永久正电荷"},
    "Sulfonium":          {"charge": +1, "note": "锍盐永久正电荷"},
    "Pyridinium":        {"charge": +1, "note": "吡啶𬭩永久正电荷"},
    "TMPP":              {"charge": +1, "note": "TMPP 膦𬭩试剂标签，固定 +1"},
    "Phosphonium":       {"charge": +1, "note": "膦𬭩固定正电荷"},
}


class IonError(ValueError):
    """离子化参数不合法。"""


@dataclass(frozen=True)
class IonSpecies:
    """一个具体离子：中性质量 + 总电荷 + 加合模式 + 固定电荷 + 质子数。"""

    neutral_mass: float
    charge: int               # 总电荷（带符号）
    adduct: str = "H"
    fixed_charge: int = 0     # 其中固定电荷部分（带符号）

    def __post_init__(self) -> None:
        if self.charge == 0:
            raise IonError("电荷不能为 0；中性分子请用 neutral_mass 直接比较")
        if self.adduct not in ADDUCTS:
            raise IonError(f"未登记加合 {self.adduct!r}；可选 {sorted(ADDUCTS)}")
        if abs(self.fixed_charge) > abs(self.charge):
            raise IonError(f"固定电荷 {self.fixed_charge:+d} 超过总电荷 {self.charge:+d}")
        if self.fixed_charge and self.adduct != "H":
            raise IonError("固定电荷离子只与质子化/去质子化组合建模；加 Na/K/NH4 时的固定电荷需单独核对")

    @property
    def protons(self) -> int:
        """质子（或去质子）个数：总电荷里扣除固定电荷后的部分。"""
        return self.charge - self.fixed_charge

    def mz(self) -> float:
        """理论 m/z（恒为正数，电荷符号见 charge 字段）。

        固定电荷部分不再加质子，这正是"独立建模"的核心：
        例如固定 +1 的季铵肽 [M]⁺ 的 m/z = M/1，而不是 (M+H)/1。
        """
        z = abs(self.charge)
        if self.adduct == "H":
            ion_mass = self.neutral_mass + self.protons * PROTON
        elif self.adduct == "H-neg":
            ion_mass = self.neutral_mass + self.protons * PROTON
        else:
            ion_mass = self.neutral_mass + z * ADDUCTS[self.adduct]["mass"]
        return ion_mass / z

    def describe(self) -> str:
        ad = ADDUCTS[self.adduct]["label"]
        bits = [f"z={self.charge:+d}", ad.replace("n", str(abs(self.charge)))]
        if self.fixed_charge:
            bits.append(f"其中固定电荷 {self.fixed_charge:+d}（质子数 {self.protons:+d}）")
        return "；".join(bits)


def mz_of(neutral_mass: float, charge: int, adduct: str = "H",
          fixed_charge: int = 0) -> float:
    """便捷函数：中性质量 → 理论 m/z。"""
    return IonSpecies(neutral_mass, charge, adduct, fixed_charge).mz()


def neutral_from_ion(mz: float, charge: int, adduct: str = "H",
                     fixed_charge: int = 0) -> float:
    """m/z → 中性质量（IonSpecies.mz 的逆运算，固定电荷按同一约定处理）。"""
    sp = IonSpecies(0.0, charge, adduct, fixed_charge)
    z = abs(charge)
    ion = mz * z
    if adduct in ("H", "H-neg"):
        return ion - sp.protons * PROTON
    return ion - z * ADDUCTS[adduct]["mass"]


def intrinsic_charge(sequence: str, ph: float = 2.0,
                     cterm_amide: bool = False) -> int:
    """近似净电荷（形式电荷和，不是精确 pKa 滴定）。

    规则：某基团 pKa 与 pH 的关系决定是否带电荷（pH < pKa 的碱带 +1，pH > pKa 的酸带 −1）。
    只作"这肽在酸性条件下大概几个电荷"的提示，用于选 z 值排查，不写进报告结论。
    """
    if ph <= 0:
        raise IonError("pH 必须为正")
    q = 0
    q += 1 if ph < _PKA["Nterm"] else 0
    if not cterm_amide:
        q -= 1 if ph > _PKA["Cterm"] else 0
    for ch in sequence.upper():
        if ch == "K" and ph < _PKA["K"]:
            q += 1
        elif ch == "R" and ph < _PKA["R"]:
            q += 1
        elif ch == "H" and ph < _PKA["H"]:
            q += 1
        elif ch == "D" and ph > _PKA["D"]:
            q -= 1
        elif ch == "E" and ph > _PKA["E"]:
            q -= 1
    return q


def delta_mz_to_neutral(delta_mz: float, charge: int) -> float:
    """两个实测 m/z 之差 → 中性质量差。必须指定电荷，否则会差整数倍。"""
    if charge == 0:
        raise IonError("电荷不能为 0")
    return delta_mz * abs(charge)


def propagate_delta_error(err1_da: float, err2_da: float) -> float:
    """两个实测值各自的不确定度 → 差值的不确定度（独立误差，方和根）。"""
    return (err1_da ** 2 + err2_da ** 2) ** 0.5


def fixed_charge_groups() -> dict[str, dict]:
    return dict(FIXED_CHARGE_GROUPS)