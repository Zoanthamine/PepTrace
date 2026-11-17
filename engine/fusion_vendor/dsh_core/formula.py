"""单同位素质量与元素式代数（DSH 版独立实现）。

设计要点
- 原子量表**独立手打**（IUPAC/CIAAW 常用值），不复用项目里已有的 `mass_ledger.py`，
  这样两边算出的数如果一致，才算真正互证；不一致就是有 bug，必须查。
- `Formula` 是 `dict[str, int]` 的薄包装，支持 +、-、* 和比较；空式等于零。
- 所有质量都是**中性单同位素质量**（Da）。带电物种必须另外用 `mz()` 换算，
  绝不能把乘了电荷的质量当中性质量写进知识库（这是书里几次踩过的坑）。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

# CIAAW/IUPAC 常用单同位素原子质量（Da）。顺序与取值刻意与项目旧脚本独立。
MONO: dict[str, float] = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
    "P": 30.97376199842,
    "F": 18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I": 126.9044719,
    "Na": 22.9897692820,
    "K": 38.9637064864,
    "Li": 7.01600455,
    "Se": 79.9165218,
    "Si": 27.97692653465,
    "B": 11.00930536,
}

# 电子质量与质子质量：m/z 换算必须用质子质量（不是 H 原子质量）。
ELECTRON = 0.000548579909065
PROTON = 1.007276466621


class FormulaError(ValueError):
    """元素式非法（未知元素、非整数下标等）。"""


def _clean(d: Mapping[str, int | float]) -> dict[str, int]:
    out: dict[str, int] = {}
    for el, n in d.items():
        if el not in MONO:
            raise FormulaError(f"未知元素符号 {el!r}；可接受：{', '.join(sorted(MONO))}")
        n_int = int(round(float(n)))
        if abs(float(n) - n_int) > 1e-9:
            raise FormulaError(f"元素 {el} 的下标必须是整数，收到 {n!r}")
        if n_int:
            out[el] = n_int
    return out


@dataclass(frozen=True)
class Formula:
    """不可变元素式。构造时自动去掉下标 0 的元素，保证 `==` 与哈希可用于查重。"""

    counts: tuple[tuple[str, int], ...]

    @staticmethod
    def of(d: Mapping[str, int | float] | None = None, **kwargs: int) -> "Formula":
        m: dict[str, int | float] = {}
        if d:
            m.update(d)
        m.update(kwargs)
        return Formula(tuple(sorted(_clean(m).items())))

    def as_dict(self) -> dict[str, int]:
        return dict(self.counts)

    def __add__(self, other: "Formula | Mapping[str, int]") -> "Formula":
        o = other if isinstance(other, Formula) else Formula.of(other)
        m = self.as_dict()
        for el, n in o.counts:
            m[el] = m.get(el, 0) + n
        return Formula.of(m)

    def __sub__(self, other: "Formula | Mapping[str, int]") -> "Formula":
        o = other if isinstance(other, Formula) else Formula.of(other)
        m = self.as_dict()
        for el, n in o.counts:
            m[el] = m.get(el, 0) - n
        return Formula.of(m)

    def __mul__(self, k: int) -> "Formula":
        return Formula.of({el: n * int(k) for el, n in self.counts})

    __rmul__ = __mul__

    def __neg__(self) -> "Formula":
        return self * -1

    def __bool__(self) -> bool:
        return bool(self.counts)

    def mass(self) -> float:
        """中性单同位素质量（Da）。"""
        return sum(MONO[el] * n for el, n in self.counts)

    def dbe(self) -> float:
        """不饱和度（等价 DBE = C − H/2 + N/2 + 1），只用于提示，不做结构判断。"""
        d = self.as_dict()
        return d.get("C", 0) - d.get("H", 0) / 2 + d.get("N", 0) / 2 + 1

    def name(self) -> str:
        """Hill 顺序式（C、H 先，其余按字母），便于人读与日志比对。"""
        return formula_str(self.as_dict())

    def __str__(self) -> str:  # pragma: no cover - 仅调试
        return self.name()


def formula_str(d: Mapping[str, int]) -> str:
    """把 {C:2,H:2,O:1} 写成 'C2H2O'（下标 1 省略，负数写成 -n）。"""
    if not d:
        return "∅"
    rest = {k: int(v) for k, v in d.items() if int(v) != 0}
    if not rest:
        return "∅"
    order: list[str] = []
    for el in ("C", "H"):
        if el in rest:
            order.append(el)
    order += sorted(el for el in rest if el not in ("C", "H"))
    parts = []
    for el in order:
        n = rest[el]
        parts.append(el if n == 1 else f"{el}{n}")
    return "".join(parts)


def delta_str(d: Mapping[str, int]) -> str:
    """带符号的净变化写法，例如 'C5H9N-1O-1' → 'C+5H+9N-1O-1' 的紧凑形式。"""
    if not d or not any(int(v) for v in d.values()):
        return "0"
    parts = []
    for el in ("C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Na", "K"):
        if el in d and int(d[el]):
            n = int(d[el])
            parts.append(f"{el}{n:+d}" if n != 1 else el)
    for el in sorted(k for k in d if k not in ("C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Na", "K")):
        n = int(d[el])
        if n:
            parts.append(f"{el}{n:+d}" if n != 1 else el)
    return "".join(parts)


def mass_of(d: Mapping[str, int]) -> float:
    return Formula.of(d).mass()


def mz(neutral_mass: float, charge: int, proton_mass: float = PROTON) -> float:
    """[M+zH]^z+ 的 m/z。charge 必须 > 0（负离子请传负电荷并用 |z| 处理加合）。"""
    if charge == 0:
        raise ValueError("电荷不能为 0")
    z = abs(int(charge))
    return (neutral_mass + z * proton_mass) / z


def neutral_from_mz(m_z: float, charge: int, proton_mass: float = PROTON) -> float:
    """由 m/z 反推中性质量（假定 [M+zH]^z+，无其他加合离子）。"""
    z = abs(int(charge))
    return m_z * z - z * proton_mass


def adduct_mz(neutral_mass: float, charge: int, adduct: str = "H") -> float:
    """通用加合离子换算。adduct 支持 H / Na / K / NH4；charge 为加合个数。"""
    z = abs(int(charge))
    cation = {
        "H": PROTON,
        "Na": MONO["Na"] - ELECTRON,
        "K": MONO["K"] - ELECTRON,
        "NH4": MONO["N"] + 4 * MONO["H"] - ELECTRON,
    }
    if adduct not in cation:
        raise ValueError(f"不支持的加合离子 {adduct!r}；可用：{', '.join(cation)}")
    return (neutral_mass + z * cation[adduct]) / z


def sum_formulas(items: Iterable[Formula | Mapping[str, int]]) -> Formula:
    total = Formula.of()
    for it in items:
        total = total + (it if isinstance(it, Formula) else Formula.of(it))
    return total