"""多肽状态模型：序列 + 端基 + 修饰 → 唯一规范态 → 分子式/中性质量。

为什么要有"规范态"
------------------
同一件事可以有多种写法：`cterm="amide"` 和 `cterm="acid"` + C 端 `C-amide` 修饰
在化学上完全等价。如果不在入口把它们折叠成同一个对象，下游就会出现
"质量一样但候选不一样"的幽灵差异（这正是执行计划第 4 节第 4 条要固化的失败案例）。
所以：**所有等价写法在这里归一，下游只认 PeptideState。**

质量约定
--------
- 一切质量都是**中性单同位素质量**（Da）。
- 游离肽链（N 端游离、C 端游离酸）= Σ 残基 + H2O。
- N 端酰化/烷氧羰基化 = 在游离胺上加共价净增量（组件库里的 delta 就是"相对游离胺的净增量"）。
- C 端酰胺 = +N +H −O（把 COOH 的 OH 换成 NH2）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .components import ComponentLibrary
from .formula import Formula

# ---------------------------------------------------------------- 别名表
# 端基/修饰组的标准名 → 组件库中的正式键。用户输入常见写法都收敛到这里。
GROUP_ALIASES: dict[str, str] = {
    "acetyl": "Ac",
    "ac": "Ac",
    "ac2o": "Ac",
    "nh2": "C-amide",
    "amide": "C-amide",
    "conh2": "C-amide",
    "ome": "C-methyl-ester",
    "coome": "C-methyl-ester",
    "methyl ester": "C-methyl-ester",
    "cf3co": "TFA",
    "tfa-salt": "TFA",
    "fdn": "Formyl",
    "cho": "Formyl",
    "pbf": "Pbf",
    "pmc": "Pmc",
    "boc": "Boc",
    "fmoc": "Fmoc",
    "cbz": "Cbz",
    "z": "Cbz",
    "alloc": "Alloc",
    "dde": "Dde",
    "acm": "Acm",
    "cam": "Carbamidomethyl",
    "carbamidomethyl": "Carbamidomethyl",
    "phospho": "Phospho",
    "p": "Phospho",
    "oxidation": "Oxidation",
    "ox": "Oxidation",
    "meto": "Oxidation",
    "disulfide": "Disulfide",
    "ss": "Disulfide",
}

# 端基写法（"free"/"none"/"" 表示游离）
FREE_TOKENS = {"free", "none", "", "h", "nh", "unmodified", "未修饰", "游离"}

# 每个基团允许出现在哪些位点；用于**硬报错**（化学上不可能的搭配）。
# 只放"不可能"的，不放"少见"的——宁可给 warning 也不能误拒合法输入。
# 每个基团允许出现在哪些侧链（T04 扩充的组件约束表）。
# 只写"化学上不可能"的搭配：宁可漏报 warning，也不能误拒合法输入
# （例如 Thr(tBu) 是合法的，早期的 _CYSMOD 写法会误拒它）。
_RESIDUE_SCOPE: dict[str, set[str]] = {
    "Pbf": {"R"},
    "Pmc": {"R"},
    "Acm": {"C"},
    "StBu": {"C"},
    "Carbamidomethyl": {"C"},
    "Carboxymethyl": {"C"},
    "Sulfo": {"C", "S", "T", "Y"},      # S-磺酸 / O-硫酸酯
    "Phospho": set("STY"),
    "tBu": {"C", "S", "T", "D", "E", "Y"},   # S-/O-叔丁基
    "Bzl": {"C", "S", "T", "D", "E", "Y"},   # S-/O-苄基
    "Dde": {"K"},
    "Alloc": {"K"},
}
_OXIDIZABLE = set("MWCH")
_ESTER_OR_AMIDE = {"C-amide", "C-methyl-ester"}

# 位点在侧链上时，这些基团必须落在带侧链胺的 Lys（胺保护基的通用规则）
_AMINE_PROTECTING = {"Fmoc", "Boc", "Cbz", "Alloc", "Dde", "Trt", "Mmt", "2-ClTrt", "Ac", "Formyl", "All", "Bzl"}


class StateError(ValueError):
    """规范态构造失败（位点/基团不可能等）。"""

    def __init__(self, code: str, message: str, hint: str = "") -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.hint = hint


@dataclass(frozen=True)
class Modification:
    """单个修饰。site 归一为 'N' | 'C' | 'r<1-based 序号>'。"""

    site: str
    group: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.site}/{self.group}"


def canonical_group(name: str) -> str:
    key = name.strip()
    if key in GROUP_ALIASES:
        return GROUP_ALIASES[key]
    low = key.lower()
    for alias, canon in GROUP_ALIASES.items():
        if alias == low:
            return canon
    return key


def canonical_site(site: str | int, seq_len: int) -> str:
    """把位点写法归一为 'N' / 'C' / 'r<k>'。

    接受：'N'/'nterm'/'n-term'，'C'/'cterm'/'c-term'，整数或数字字符串（1-based），
    以及 'K3'/'R5'/'C7' 这种"残基字母+序号"（字母只用于人读，序号才是权威）。
    """
    if isinstance(site, int):
        k = site
    else:
        s = str(site).strip()
        low = s.lower()
        if low in {"n", "nterm", "n-term", "nt", "n_terminus"}:
            return "N"
        if low in {"c", "cterm", "c-term", "ct", "c_terminus"}:
            return "C"
        digits = "".join(ch for ch in s if ch.isdigit())
        if not digits:
            raise StateError("E_BAD_SITE", f"无法识别的修饰位点 {site!r}", "位点写 'N'、'C' 或 1-based 残基序号（如 3 或 'K3'）")
        k = int(digits)
    if k < 1 or k > seq_len:
        raise StateError(
            "E_SITE_OUT_OF_RANGE",
            f"修饰位点 {site!r} 超出序列长度（序列 {seq_len} 个残基）",
            "位点从 1 开始计数",
        )
    return f"r{k}"


@dataclass
class PeptideState:
    """规范态：序列 + N 端基团 + C 端是否酰胺 + 侧链修饰列表。"""

    sequence: str
    nterm_group: str | None = None
    cterm_amide: bool = False
    modifications: tuple[Modification, ...] = ()
    # 输入里出现过但被折叠掉的等价写法，用于报告里说明"我按什么理解你的输入"
    folded: tuple[str, ...] = field(default_factory=tuple)

    # ---------------- 序列信息 ----------------
    @property
    def length(self) -> int:
        return len(self.sequence)

    def residues(self) -> list[str]:
        return list(self.sequence)

    def has(self, aa: str) -> bool:
        return aa.upper() in self.sequence.upper()

    def positions_of(self, aa: str) -> list[int]:
        aa = aa.upper()
        return [i + 1 for i, ch in enumerate(self.sequence.upper()) if ch == aa]

    # ---------------- 质量 ----------------
    def formula(self, lib: ComponentLibrary) -> Formula:
        total = Formula.of({"H": 2, "O": 1})  # 游离肽 + H2O
        for aa in self.sequence:
            total = total + lib.residue(aa.upper()).formula
        if self.nterm_group:
            total = total + lib.delta(self.nterm_group).formula
        for mod in self.modifications:
            total = total + lib.delta(mod.group).formula
        if self.cterm_amide:
            # COOH → CONH2：+N +H −O
            total = total + Formula.of({"N": 1, "H": 1, "O": -1})
        return total

    def mass(self, lib: ComponentLibrary) -> float:
        return self.formula(lib).mass()

    def describe(self) -> str:
        """一句话描述规范态，写进报告，避免"我以为的输入"和"实际解析的输入"不一致。"""
        bits = [self.sequence]
        if self.nterm_group:
            bits.append(f"N端={self.nterm_group}")
        else:
            bits.append("N端游离")
        bits.append("C端酰胺" if self.cterm_amide else "C端游离酸")
        if self.modifications:
            bits.append("修饰=" + ",".join(str(m) for m in self.modifications))
        return " | ".join(bits)


def build_state(
    sequence: str,
    nterm: str | None,
    cterm: str | None,
    modifications: Iterable[dict | Modification],
    lib: ComponentLibrary,
) -> PeptideState:
    """构造规范态，并把等价写法折叠掉（C 端酰胺、N 端酰基）。"""
    seq = "".join(ch for ch in str(sequence).strip() if not ch.isspace() and ch not in "-_")
    if not seq:
        raise StateError("E_BAD_SEQUENCE", "序列为空", "给出一字母序列，如 'GC' 或 'Ac-WG-NH2'")
    bad = [ch for ch in seq if ch.upper() not in lib.residues]
    if bad:
        raise StateError(
            "E_BAD_SEQUENCE",
            f"序列含未登记残基 {''.join(bad)!r}",
            f"已登记残基：{''.join(sorted(lib.residues))}；三字母写法请先转成一字母",
        )

    folded: list[str] = []

    # --- N 端：'free' / 组名 ---
    nterm_group: str | None = None
    if nterm is not None and str(nterm).strip():
        tok = str(nterm).strip()
        if tok.lower() in FREE_TOKENS:
            nterm_group = None
        else:
            g = canonical_group(tok)
            if not lib.has_delta(g):
                raise StateError("E_UNKNOWN_GROUP", f"未登记 N 端基团 {tok!r}", "检查拼写，或在 components 里新增该基团")
            nterm_group = g

    # --- C 端：'acid' / 'amide' ---
    cterm_amide = False
    if cterm is not None and str(cterm).strip():
        tok = str(cterm).strip().lower()
        if tok in {"acid", "free", "cooh", "oh", "-oh", "酸", "游离"}:
            cterm_amide = False
        elif tok in {"amide", "nh2", "conh2", "酰胺", "amidated"}:
            cterm_amide = True
        else:
            raise StateError("E_BAD_CTERM", f"C 端写法 {cterm!r} 无法识别", "用 'acid' 或 'amide'")

    # --- 修饰 ---
    mods: list[Modification] = []
    for raw in modifications:
        if isinstance(raw, Modification):
            mods.append(raw)
            continue
        if "group" not in raw:
            raise StateError("E_BAD_MODIFICATION", f"修饰缺少 group 字段：{raw!r}", '形如 {"site": "N", "group": "Fmoc"}')
        g = canonical_group(str(raw["group"]))
        if not lib.has_delta(g):
            raise StateError("E_UNKNOWN_GROUP", f"未登记修饰基团 {raw['group']!r}", "检查拼写；拼错和未知字段都会被报错，不会静默忽略")
        site = canonical_site(raw.get("site", "N"), len(seq))
        if g in _ESTER_OR_AMIDE:
            if site != "C":
                raise StateError("E_INCOMPATIBLE_SITE", f"基团 {g} 只能加在 C 端，收到位点 {site}", "C 端修饰请写 site='C'，或直接用 cterm='amide'")
            if g == "C-amide":
                cterm_amide = True
                folded.append("C 端 C-amide 修饰 ≡ cterm=amide")
                continue
        if site in {"N", "C"} and g in _AMINE_PROTECTING and site == "C":
            raise StateError("E_INCOMPATIBLE_SITE", f"N 端保护基 {g} 不能加在 C 端", "N 端基团用 nterm 字段或 site='N'")
        if site.startswith("r"):
            k = int(site[1:])
            res = seq[k - 1].upper()
            scope = _RESIDUE_SCOPE.get(g)
            if scope and res not in scope:
                raise StateError(
                    "E_INCOMPATIBLE_SITE",
                    f"{g} 不能加在 {res} 侧链（位点 {k}）",
                    f"{g} 只适用于 {'/'.join(sorted(scope))}；确认序号是否 1-based、残基是否写错",
                )
            if g == "Oxidation" and res not in _OXIDIZABLE:
                raise StateError("E_INCOMPATIBLE_SITE", f"氧化位点 {k} 是易感残基以外的 {res}", "Met/Trp/His/Cys 才按氧化处理")
            if g in _AMINE_PROTECTING:
                # Fmoc 等也可上 Lys ε-胺，但 residue 必须是 K
                if res != "K":
                    raise StateError(
                        "E_INCOMPATIBLE_SITE",
                        f"{g} 是胺保护基，位点 {k} 是 {res}（无侧链伯胺）",
                        "只能加在 N 端或 Lys 侧链",
                    )
        mods.append(Modification(site, g))

    if nterm_group and nterm_group in _ESTER_OR_AMIDE:
        raise StateError("E_INCOMPATIBLE_SITE", f"{nterm_group} 不能作为 N 端基团", "C 端修饰请用 cterm 字段")

    # 统一折叠：site='N' 的修饰就是 N 端封端，归入 nterm_group（质量相同，但规范态必须唯一）
    n_site = [m for m in mods if m.site == "N"]
    if len(n_site) > 1:
        raise StateError(
            "E_DUPLICATE_MOD",
            f"N 端修饰重复出现 {len(n_site)} 次（{', '.join(m.group for m in n_site)}）",
            "同一基团一次性加上即可；重复保护不会叠加质量",
        )
    if n_site:
        if len({m.group for m in n_site}) > 1:
            raise StateError(
                "E_DUPLICATE_NTERM",
                f"N 端出现多个不同基团：{', '.join(m.group for m in n_site)}",
                "N 端只能有一个封端基团",
            )
        if nterm_group:
            raise StateError(
                "E_DUPLICATE_NTERM",
                f"N 端基团重复：nterm={nterm_group}，modifications 里又有 N/{n_site[0].group}",
                "二者只写一处；推荐统一用 nterm，或统一用 modifications",
            )
        nterm_group = n_site[0].group
        mods = [m for m in mods if m.site != "N"]
        folded.append(f"N 端修饰 {nterm_group} ≡ nterm={nterm_group}")

    # 去重：同一 (site, group) 出现两次是明确的输入错误
    seen: dict[tuple[str, str], int] = {}
    for m in mods:
        key = (m.site, m.group)
        seen[key] = seen.get(key, 0) + 1
    dup = [f"{k[0]}/{k[1]}×{v}" for k, v in seen.items() if v > 1]
    if dup:
        raise StateError("E_DUPLICATE_MOD", f"同一修饰重复出现：{', '.join(dup)}", "同一基团一次性加上即可；若确实要加两次（如两个 Lys 各加一个），请分别写不同位点")

    return PeptideState(
        sequence=seq.upper(),
        nterm_group=nterm_group,
        cterm_amide=cterm_amide,
        modifications=tuple(sorted(mods, key=lambda m: (m.site, m.group))),
        folded=tuple(folded),
    )