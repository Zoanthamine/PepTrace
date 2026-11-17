"""同位素与平均质量模型（T04）。

为什么需要
----------
执行计划第 6 节要求：质量组件必须区分**名义 / 平均 / 单同位素**，低分辨检索要按实际误差区间，
不能"四舍五入当平均质量"。只有单同位素质量时，低分辨数据（尤其 MALDI、四极杆低分辨）会系统偏小；
反过来，用平均质量去比高分辨单同位素峰又是错的。所以这里同时提供：
  * 单同位素质量（formula.Formula.mass()，另一套原子量表）
  * 平均质量（本模块，IUPAC 平均原子量）
  * 同位素包络 M、M+1、M+2（本模块，按元素天然丰度卷积）

数据来源与精度
--------------
平均原子量取 IUPAC 常用值（H 1.00794、C 12.0107、N 14.0067、O 15.9994、S 32.065、
P 30.973762 等）。同位素丰度取 NIST 常用值。这是**近似模型**：用于低分辨/平均质量口径与
同位素峰位提示，不用于高分辨精确质量的替代（高分辨一律用单同位素 + 净式）。
"""
from __future__ import annotations

from .formula import Formula

# ---------------------------------------------------------------- 平均原子量
AVERAGE: dict[str, float] = {
    "H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "S": 32.065,
    "P": 30.973762, "F": 18.9984032, "Cl": 35.453, "Br": 79.904, "I": 126.90447,
    "Na": 22.98976928, "K": 39.0983, "Si": 28.0855, "Se": 78.96, "B": 10.811,
    "Li": 6.941, "Fe": 55.845, "Zn": 65.38,
}

# ---------------------------------------------------------------- 同位素模式
# 每项：(相对单同位素的质量偏移 Da, 天然丰度)。只保留到 M+2 的主峰。
ISOTOPES: dict[str, tuple[tuple[float, float], ...]] = {
    "H": ((0.0, 0.999885), (1.007825032, 0.000115)),
    "C": ((0.0, 0.9893), (1.003354838, 0.0107)),
    "N": ((0.0, 0.99632), (0.997034893, 0.00368)),
    "O": ((0.0, 0.99757), (1.004217136, 0.00038), (2.004244992, 0.00205)),
    "S": ((0.0, 0.9499), (1.9957959, 0.0425), (0.999387, 0.0075), (3.995799, 0.0001)),
    "P": ((0.0, 1.0),),
    "Cl": ((0.0, 0.7576), (1.99704992, 0.2424)),
    "Br": ((0.0, 0.5069), (1.9979535, 0.4931)),
    "I": ((0.0, 1.0),),
    "F": ((0.0, 1.0),),
    "Na": ((0.0, 1.0),),
    "K": ((0.0, 0.932581), (1.9980111, 0.000117), (1.9980111, 0.067302)),
    "Se": ((0.0, 0.4961), (1.9978, 0.0923), (3.9900, 0.2377), (1.9948, 0.0763)),
}


def average_mass_of(formula: Formula) -> float:
    """平均质量（Da）。未知元素抛 KeyError，不静默按 0 处理。"""
    total = 0.0
    unknown = []
    for el, n in formula.counts:
        if el not in AVERAGE:
            unknown.append(el)
            continue
        total += AVERAGE[el] * n
    if unknown:
        raise KeyError(f"平均原子量表中没有元素 {sorted(set(unknown))}，不能给平均质量")
    return total


def average_mass_of_counts(counts: dict[str, int]) -> float:
    return average_mass_of(Formula.of(counts))


def _pow_pattern(pattern: tuple[tuple[float, float], ...], n: int,
                 max_peaks: int = 4) -> list[tuple[float, float]]:
    """n 个该元素原子的同位素分布（截断到 max_peaks 个峰）。"""
    if n <= 0:
        return [(0.0, 1.0)]
    dist = [(0.0, 1.0)]
    for _ in range(n):
        new: dict[int, list[float]] = {}
        for m1, p1 in dist:
            for m2, p2 in pattern:
                key = int(round(m1 + m2))          # 按名义增量分箱
                slot = new.setdefault(key, [0.0, 0.0])
                slot[0] += p1 * p2                 # 概率
                slot[1] += p1 * p2 * (m1 + m2)     # 质量加权，用于算箱内平均偏移
        dist = [(slot[1] / slot[0], slot[0]) for key, slot in sorted(new.items())
                if slot[0] > 0]
        dist = [(m, p) for m, p in dist if p > 1e-12]
        dist.sort(key=lambda x: -x[1])
        dist = sorted(dist[:max_peaks])
    return dist


def envelope(formula: Formula, max_offset: int = 2) -> list[dict]:
    """同位素包络：返回 [{offset, mass_shift, probability, rel_intensity}, ...]。

    offset=0 即单同位素峰。概率是绝对丰度（同一分子各峰概率和 ≲ 1，截断后略小）。
    """
    dist: list[tuple[float, float]] = [(0.0, 1.0)]
    for el, n in formula.counts:
        pat = ISOTOPES.get(el)
        if pat is None:
            raise KeyError(f"同位素表中没有元素 {el}")
        if n == 0:
            continue
        el_dist = _pow_pattern(pat, n)
        new: dict[int, list[float]] = {}
        for m1, p1 in dist:
            for m2, p2 in el_dist:
                key = int(round(m1 + m2))
                slot = new.setdefault(key, [0.0, 0.0])
                slot[0] += p1 * p2
                slot[1] += p1 * p2 * (m1 + m2)
        dist = [(slot[1] / slot[0], slot[0]) for key, slot in sorted(new.items()) if slot[0] > 1e-12]
    dist.sort(key=lambda x: x[0])
    base = dist[0][1] if dist else 0.0
    out = []
    for shift, prob in dist:
        offset = int(round(shift))
        if offset > max_offset:
            continue
        out.append({
            "offset": offset,
            "mass_shift": round(shift, 6),
            "probability": round(prob, 8),
            "rel_intensity": round(prob / base, 6) if base else None,
        })
    return out


def monoisotopic(formula: Formula) -> float:
    return formula.mass()