"""质量匹配：只做数值判断，不掺条件判断（条件在 candidates.py）。

三个口径，互不混淆
------------------
1. `delta`：实测/理论**中性质量差**（产物 − 起始物），单位 Da。
2. `neutral_mass`：实测中性单同位素质量。
3. `mz`：实测 m/z，必须先按电荷与加合离子换成中性质量再比。

低分辨（只有整数 amu 的规则）单独成档：整数规则只能匹配到"±0.5 amu 的箱"，
绝不能把它算成 δ=0 的精确匹配——那是把低分辨信息伪装成高分辨结论。
"""
from __future__ import annotations

from dataclasses import dataclass

from .rules import Rule


@dataclass
class MassMatch:
    rule: Rule
    predicted: float          # 预测的中性质量差（Da）
    observed: float           # 观测（已归一到中性质量差）
    error_da: float           # 实测 − 理论（全文统一这个方向）
    error_ppm: float | None
    grade: str                # exact | nominal | nominal_multiple
    unit_ratio: float = 1.0   # nominal 倍数情形下观测/单位 的比值

    @property
    def abs_error(self) -> float:
        return abs(self.error_da)


def ppm_error(observed: float, predicted: float) -> float | None:
    if predicted == 0:
        return None
    return (observed - predicted) / abs(predicted) * 1e6


def match_delta(rules: list[Rule], observed_delta: float, tol_da: float) -> list[MassMatch]:
    """按质量差匹配规则。tol_da 是容差（Da）；整数规则另加 0.5 amu 的箱宽。

    "及其倍数"的规则（如 SO3 逐级磺化 +80/+160/+240）单独成档：
    整数规则标 `nominal_multiple`，精确质量规则标 `exact_multiple`，绝不当成单次加成的精确匹配。
    """
    out: list[MassMatch] = []
    for r in rules:
        # ---- 可重复单元 ----
        if r.nominal_multiples and (r.delta_mass is not None or r.nominal is not None):
            base = r.delta_mass if r.delta_mass is not None else r.nominal
            if not base:
                continue
            exact = r.delta_mass is not None
            best: MassMatch | None = None
            for k in range(1, 9):
                pred = base * k
                window = tol_da * k if exact else (0.5 * k + tol_da * k)
                err = observed_delta - pred
                if abs(err) <= window:
                    grade = "exact_multiple" if exact else "nominal_multiple"
                    cand = MassMatch(r, pred, observed_delta, err, ppm_error(observed_delta, pred), grade, k)
                    if best is None or cand.abs_error < best.abs_error:
                        best = cand
            if best:
                out.append(best)
            continue
        # ---- 仅整数 amu ----
        if r.is_low_res_only and r.nominal is not None:
            err = observed_delta - r.nominal
            if abs(err) <= 0.5 + tol_da:
                out.append(MassMatch(r, r.nominal, observed_delta, err, ppm_error(observed_delta, r.nominal), "nominal"))
            continue
        # ---- 精确质量差 ----
        if r.delta_mass is None:
            continue
        err = observed_delta - r.delta_mass
        if abs(err) <= tol_da:
            out.append(MassMatch(r, r.delta_mass, observed_delta, err, ppm_error(observed_delta, r.delta_mass), "exact"))
    return out


def match_absolute(rules: list[Rule], target_mass: float, observed_neutral: float, tol_da: float) -> list[MassMatch]:
    """按绝对中性质量匹配：预测产物质量 = 目标质量 + 规则质量差。"""
    implied = observed_neutral - target_mass
    out: list[MassMatch] = []
    for r in rules:
        if r.delta_mass is None and r.nominal is None:
            continue
        base = r.delta_mass if r.delta_mass is not None else r.nominal
        assert base is not None
        pred_product = target_mass + base
        err = observed_neutral - pred_product
        grade = "exact" if r.delta_mass is not None else "nominal"
        window = tol_da if r.delta_mass is not None else 0.5 + tol_da
        if abs(err) <= window:
            out.append(MassMatch(r, base, implied, err, ppm_error(observed_neutral, pred_product), grade))
    return out


def _group_key(rule) -> str:
    """同式分组键：优先净式，其次整数/倍数口径，最后按数值质量差。

    注意「异构不变质量」类规则（28 条）净式与 nominal 都是 None、只有 delta_mass=0，
    早期实现在这里会 `None.__format__` 崩掉；分组键必须把它们也归好。
    """
    if rule.delta is not None:
        return rule.delta.name()
    if rule.nominal is not None:
        return f"nominal{rule.nominal:+.0f}"
    if rule.delta_mass is not None:
        return f"ΔM={rule.delta_mass:+.6f}"
    return "无可计算质量差"


def group_same_delta(matches: list[MassMatch]) -> list[dict]:
    """同式/同质量候选分组：净式相同但机理、位点不同的规则必须并列显示。"""
    groups: dict[str, dict] = {}
    for m in matches:
        key = _group_key(m.rule)
        g = groups.setdefault(key, {"delta_formula": key, "delta_mass": m.rule.delta_mass, "rule_ids": [], "names": []})
        g["rule_ids"].append(m.rule.rule_id)
        g["names"].append(m.rule.name_cn)
    out = [g for g in groups.values() if len(g["rule_ids"]) > 1]
    out.sort(key=lambda g: -len(g["rule_ids"]))
    return out