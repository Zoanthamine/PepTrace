"""T07 输入/输出适配器：把外部数据（LC–MS 前体、峰对、MS/MS 峰表）转成病例输入，
以及把结果 JSON 摊平成表格。

设计约束：
- 适配器只做「形状转换」，不猜化学。任何缺字段都写进 note 而不是默认成看似有依据的值。
- 转换结果仍是普通 dict，必须再过 dsh_core.schema.parse_case 才算合法输入，
  这样所有校验只有一处实现。
- 二次谱（b/y）不参与候选打分，只作验证证据，见 dsh_core/msms_assist.py（若可用）。
"""

from __future__ import annotations

from typing import Any, Iterable

__all__ = [
    "case_from_precursor",
    "case_from_peak_pair",
    "case_from_neutral_mass",
    "case_from_msms_scan",
    "rows_from_result",
    "tsv_from_result",
]


def _clean(d: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in d.items() if v is not None}


def case_from_precursor(
    *,
    sequence: str,
    precursor_mz: float,
    charge: int,
    adduct: str = "H",
    tol_da: float | None = None,
    tol_ppm: float | None = None,
    case_id: str | None = None,
    fixed_charge: int = 0,
    nterm: str | None = None,
    cterm: str | None = None,
    modifications: list[dict[str, Any]] | None = None,
    conditions: dict[str, Any] | None = None,
    notes: str | None = None,
) -> dict[str, Any]:
    """从 LC–MS 前体离子（m/z + 电荷 + 加合物）构造病例。

    m/z 会被引擎按「加合物 + 电荷」换回中性质量再与目标质量比较，
    因此调用者不需要自己减质子；把 z 与 adduct 传对即可。
    """
    observed: dict[str, Any] = {
        "kind": "mz",
        "value": float(precursor_mz),
        "z": int(charge),
        "adduct": adduct,
    }
    if tol_da is not None:
        observed["tol_da"] = float(tol_da)
    if tol_ppm is not None:
        observed["tol_ppm"] = float(tol_ppm)
    if notes:
        observed["note"] = notes
    return _clean(
        {
            "case_id": case_id,
            "sequence": sequence,
            "nterm": nterm,
            "cterm": cterm,
            "modifications": modifications,
            "fixed_charge": fixed_charge or None,
            "conditions": conditions,
            "observed": observed,
        }
    )


def case_from_peak_pair(
    *,
    sequence: str,
    mz1: float,
    mz2: float,
    z1: int = 1,
    z2: int = 1,
    adduct1: str = "H",
    adduct2: str = "H",
    tol_da: float | None = None,
    tol_ppm: float | None = None,
    case_id: str | None = None,
    fixed_charge: int = 0,
    nterm: str | None = None,
    cterm: str | None = None,
    modifications: list[dict[str, Any]] | None = None,
    conditions: dict[str, Any] | None = None,
    notes: str | None = None,
) -> dict[str, Any]:
    """从两个实测峰构造「同一分子的两个状态」的质量差病例。

    误差按两峰 ppm 方和根传播（由 schema.parse_observed 实现），
    电荷与加合物分别核对后才知道两个峰是不是同一个分子——这一步由人负责，
    适配器只把 z/adduct 原样带上，不做电荷推测。
    """
    observed: dict[str, Any] = {
        "kind": "mz_pair",
        "value": float(mz1),
        "z": int(z1),
        "adduct": adduct1,
        "value2": float(mz2),
        "z2": int(z2),
    }
    if adduct2 != adduct1:
        observed["adduct2"] = adduct2
    if tol_da is not None:
        observed["tol_da"] = float(tol_da)
    if tol_ppm is not None:
        observed["tol_ppm"] = float(tol_ppm)
    if notes:
        observed["note"] = notes
    return _clean(
        {
            "case_id": case_id,
            "sequence": sequence,
            "nterm": nterm,
            "cterm": cterm,
            "modifications": modifications,
            "fixed_charge": fixed_charge or None,
            "conditions": conditions,
            "observed": observed,
        }
    )


def case_from_neutral_mass(
    *,
    sequence: str,
    observed_neutral_mass: float,
    tol_da: float | None = None,
    tol_ppm: float | None = None,
    case_id: str | None = None,
    fixed_charge: int = 0,
    nterm: str | None = None,
    cterm: str | None = None,
    modifications: list[dict[str, Any]] | None = None,
    conditions: dict[str, Any] | None = None,
    notes: str | None = None,
) -> dict[str, Any]:
    """从 MALDI/低分辨测得的中性质量构造绝对质量病例（观测口径 = 实测 − 目标）。"""
    observed: dict[str, Any] = {"kind": "neutral_mass", "value": float(observed_neutral_mass)}
    if tol_da is not None:
        observed["tol_da"] = float(tol_da)
    if tol_ppm is not None:
        observed["tol_ppm"] = float(tol_ppm)
    if notes:
        observed["note"] = notes
    return _clean(
        {
            "case_id": case_id,
            "sequence": sequence,
            "nterm": nterm,
            "cterm": cterm,
            "modifications": modifications,
            "fixed_charge": fixed_charge or None,
            "conditions": conditions,
            "observed": observed,
        }
    )


def case_from_msms_scan(
    *,
    sequence: str,
    precursor_mz: float,
    precursor_charge: int,
    peaks: Iterable[tuple[float, float]],
    adduct: str = "H",
    tol_ppm: float = 20.0,
    case_id: str | None = None,
    nterm: str | None = None,
    cterm: str | None = None,
    modifications: list[dict[str, Any]] | None = None,
    conditions: dict[str, Any] | None = None,
    annotate: bool = True,
) -> dict[str, Any]:
    """从一次 MS/MS 扫描构造病例，并（可选）附带 b/y 归属证据。

    返回 dict 形如 ``{"case": <病例 dict>, "evidence": {...}|None}``。
    b/y 归属只用于验证候选，不参与质量匹配打分。
    """
    case = case_from_precursor(
        sequence=sequence,
        precursor_mz=precursor_mz,
        charge=precursor_charge,
        adduct=adduct,
        tol_ppm=tol_ppm,
        case_id=case_id,
        nterm=nterm,
        cterm=cterm,
        modifications=modifications,
        conditions=conditions,
        notes="由 MS/MS 扫描适配；碎片归属仅作验证证据",
    )
    evidence: dict[str, Any] | None = None
    if annotate:
        try:
            import sys
            from pathlib import Path

            tools_dir = Path(__file__).resolve().parent.parent / "tools"
            if str(tools_dir) not in sys.path:
                sys.path.insert(0, str(tools_dir))
            import msms_assist  # 依赖 pyOpenMS，缺失则降级

            evidence = {
                "available": True,
                "annotations": msms_assist.annotate_sequence(
                    sequence,
                    list(peaks),
                    precursor_charge,
                    tol_ppm,
                ),
                "note": "b/y 归属仅作候选验证证据，不参与质量匹配打分",
            }
        except Exception as exc:  # noqa: BLE001 - 缺少 pyOpenMS 时只降级，不影响主流程
            evidence = {"available": False, "reason": f"b/y 归属不可用：{type(exc).__name__}: {exc}"}
    return {"case": case, "evidence": evidence}


def rows_from_result(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """把结果 JSON 摊平成一张表：单事件候选 + 组合候选，行内用 kind 区分。"""
    rows: list[dict[str, Any]] = []
    for c in payload.get("candidates") or []:
        rows.append(
            {
                "kind": "single",
                "verdict": c.get("verdict"),
                "rule_or_key": c.get("rule_id"),
                "name": c.get("name_cn"),
                "predicted_delta_da": c.get("predicted_delta_da"),
                "error_da": c.get("error_da"),
                "error_ppm": c.get("error_ppm"),
                "mass_grade": c.get("mass_grade"),
                "condition_state": c.get("condition_state"),
                "sites": ";".join(str(s) for s in (c.get("sites") or [])),
                "excluded_reason": c.get("excluded_reason") or "",
            }
        )
    for c in payload.get("excluded") or []:
        rows.append(
            {
                "kind": "single_excluded",
                "verdict": c.get("verdict"),
                "rule_or_key": c.get("rule_id"),
                "name": c.get("name_cn"),
                "predicted_delta_da": c.get("predicted_delta_da"),
                "error_da": c.get("error_da"),
                "error_ppm": c.get("error_ppm"),
                "mass_grade": c.get("mass_grade"),
                "condition_state": c.get("condition_state"),
                "sites": ";".join(str(s) for s in (c.get("sites") or [])),
                "excluded_reason": c.get("excluded_reason") or "",
            }
        )
    for c in payload.get("combinations") or []:
        rows.append(
            {
                "kind": "combination",
                "verdict": c.get("verdict"),
                "rule_or_key": c.get("combination_key"),
                "name": " + ".join(
                    str(e.get("name_cn") or e.get("rule_id")) for e in (c.get("events") or [])
                ),
                "predicted_delta_da": c.get("predicted_delta_da"),
                "error_da": c.get("error_da"),
                "error_ppm": c.get("error_ppm"),
                "mass_grade": c.get("mass_grade"),
                "condition_state": c.get("condition_state"),
                "sites": ";".join(str(s) for s in (c.get("sites") or [])),
                "excluded_reason": c.get("excluded_reason") or "",
            }
        )
    return rows


_COLUMNS = [
    "kind",
    "verdict",
    "rule_or_key",
    "name",
    "predicted_delta_da",
    "error_da",
    "error_ppm",
    "mass_grade",
    "condition_state",
    "sites",
    "excluded_reason",
]


def tsv_from_result(payload: dict[str, Any]) -> str:
    """把结果摊平成 TSV 文本（Excel 可直接打开），列固定为 _COLUMNS。"""
    lines = ["\t".join(_COLUMNS)]
    for row in rows_from_result(payload):
        lines.append("\t".join("" if row[c] is None else str(row[c]) for c in _COLUMNS))
    return "\n".join(lines) + "\n"