"""输入契约：解析 + 严格校验。错误必须"说人话、给出处、给改法"。

设计原则（对应执行计划第 4 节第 5 条）
--------------------------------------
输入字段**严格白名单**：拼错的键（比如 `modification` 写成单数、`grup` 拼错）
一律抛 `InputError`，并给出"你是不是想写 X"。静默忽略未知字段是最危险的失败模式：
用户以为提示生效了，系统却把线索扔了。
"""
from __future__ import annotations

import difflib
import re
from dataclasses import dataclass, field
from typing import Any

from .components import ComponentLibrary
from .state import PeptideState, StateError, canonical_group, build_state


class InputError(ValueError):
    """输入不合法。带 code / field / hint，方便测试断言与报告呈现。"""

    def __init__(self, code: str, field_name: str, message: str, hint: str = "") -> None:
        super().__init__(f"[{code}] {field_name}: {message}")
        self.code = code
        self.field = field_name
        self.message = message
        self.hint = hint

    def as_dict(self) -> dict:
        return {"code": self.code, "field": self.field, "message": self.message, "hint": self.hint}


# ---------------------------------------------------------------- 试剂同义词
# 键：归一化后的小写写法；值：标准试剂名。
REAGENT_ALIASES: dict[str, str] = {
    "edt": "EDT",
    "1,2-ethanedithiol": "EDT",
    "ethanedithiol": "EDT",
    "tfa": "TFA",
    "trifluoroacetic acid": "TFA",
    "piperidine": "Piperidine",
    "pip": "Piperidine",
    "4-methylpiperidine": "4-Methylpiperidine",
    "4mp": "4-Methylpiperidine",
    "4-methylpiperidine": "4-Methylpiperidine",
    "tis": "TIS",
    "triisopropylsilane": "TIS",
    "dtt": "DTT",
    "dithiothreitol": "DTT",
    "tcep": "TCEP",
    "hf": "HF",
    "tfmsa": "TFMSA",
    "triflic acid": "TFMSA",
    "tfmsa/tfa": "TFMSA",
    "hobt": "HOBt",
    "hoat": "HOAt",
    "dic": "DIC",
    "dcc": "DCC",
    "hatu": "HATU",
    "hbtu": "HBTU",
    "pybop": "PyBOP",
    "dipea": "DIEA",
    "diea": "DIEA",
    "dipEA": "DIEA",
    "nmm": "NMM",
    "n-methylmorpholine": "NMM",
    "i2": "I2",
    "iodine": "I2",
    "dmso": "DMSO",
    "acn": "ACN",
    "acetonitrile": "ACN",
    "dcm": "DCM",
    "dichloromethane": "DCM",
    "nmp": "NMP",
    "dmf": "DMF",
    "mtbe": "MTBE",
    "ether": "Et2O",
    "diethyl ether": "Et2O",
    "methanol": "Methanol",
    "meoh": "Methanol",
    "ammonia": "Ammonia",
    "nh3": "Ammonia",
    "hydrazine": "Hydrazine",
    "nh2nh2": "Hydrazine",
    "tbaf": "TBAF",
    "zncl2": "ZnCl2",
    "znbr2": "ZnBr2",
    "pd(pph3)4": "Pd(PPh3)4",
    "water": "Water",
    "h2o": "Water",
    "phenol": "Phenol",
    "thioanisole": "Thioanisole",
    "acetone": "Acetone",
    "formaldehyde": "Formaldehyde",
    "formic acid": "FormicAcid",
    "formicacid": "FormicAcid",
    "甲酸": "FormicAcid",
    "hcooh": "FormicAcid",
    "trifluoroacetic acid": "TFA",
    "dtt": "DTT",
    "二硫苏糖醇": "DTT",
    "guanidine": "Guanidine",
    "guanidinium": "Guanidine",
    "acetic acid": "AceticAcid",
    "acoh": "AceticAcid",
    "edta": "EDTA",
    "dmap": "DMAP",
    "dmap?": "DMAP",
    "succinimide": "Succinimide",
    "nhs": "NHS",
    "hydroxylamine": "Hydroxylamine",
    "ninhydrin": "Ninhydrin",
}

# 输入里可能把"修饰"当 tag 报（例如 tags=['oxidation']）——这些不是试剂，
# 不能拿来满足任何试剂条件，只能算"用户报告了一个修饰"。
TAG_IS_MODIFICATION = {
    "oxidation", "ox", "met(o)", "deamidation", "phospho", "p", "acm",
    "carbamidomethyl", "disulfide", "ss", "cyclization", "dkp", "pyroglutamate",
    "pglu", "succinimide", "aspartimide", "isoaspartate", "isoasp",
}

TOP_FIELDS = {"case_id", "sequence", "nterm", "cterm", "modifications", "conditions",
              "observed", "notes", "fixed_charge"}
COND_FIELDS = {"reported_reagents", "denied_reagents", "tags", "stage", "notes"}
OBS_FIELDS = {"kind", "value", "z", "z2", "adduct", "tol_da", "tol_ppm", "value2", "note"}
MOD_FIELDS = {"site", "group"}


def canon_reagent(name: str) -> str | None:
    """把用户写的试剂名归一为标准名；不认识返回 None。"""
    key = re.sub(r"\s+", " ", str(name).strip().lower())
    if key in REAGENT_ALIASES:
        return REAGENT_ALIASES[key]
    # 允许 "TFA/EDT/TIS" 这种混写：逐个拆开再归一一遍
    return None


def _keep_unknown_reagent(piece: str) -> str:
    """未登记的试剂名：保留原文（首字母大写归一空白）而不是丢掉。

    政策：用户给的条件是他提供的证据。引擎不认识这个试剂，只能说"我不认识"，
    不能把它当成"用户什么也没说"——那样阴性证据会被静默忽略。
    """
    text = re.sub(r"\s+", " ", str(piece).strip())
    if not text:
        return ""
    return text[0].upper() + text[1:]


def split_reagent_blob(blob: str) -> list[str]:
    """把 'TFA/EDT/TIS(95:4:1)' 这类串拆成 ['TFA','EDT','TIS']。

    两个已修的坑：
    - 'Pd(PPh3)4' 这类含括号的催化剂名不能被括号预处理吞成 'Pd'，所以整体可识别的先整体返回；
    - 没有分隔符的串原样返回，不做括号剥离。
    """
    text = str(blob).strip()
    if not text:
        return []
    if canon_reagent(text) or canon_reagent(re.sub(r"\([^)]*\)", " ", text)):
        cleaned = re.sub(r"\([^)]*\)", " ", text) if canon_reagent(text) is None else text
        return [cleaned.strip()]
    cleaned = re.sub(r"\([^)]*\)", " ", text)
    parts = re.split(r"[/、,;，；+]|\s{2,}", cleaned)
    out = []
    for p in parts:
        p = p.strip().strip(".")
        if p:
            out.append(p)
    return out


def _suggest(key: str, allowed: set[str]) -> str:
    close = difflib.get_close_matches(key, sorted(allowed), n=1, cutoff=0.6)
    if close:
        return f"你是否想写 {close[0]!r}？"
    return f"可用字段：{', '.join(sorted(allowed))}"


def _check_fields(obj: dict, allowed: set[str], where: str) -> None:
    if not isinstance(obj, dict):
        raise InputError("E_BAD_TYPE", where, f"应为对象（dict），收到 {type(obj).__name__}")
    unknown = [k for k in obj if k not in allowed]
    if unknown:
        k = unknown[0]
        raise InputError(
            "E_UNKNOWN_FIELD",
            f"{where}.{k}",
            f"未知字段 {k!r}",
            _suggest(k, allowed),
        )


@dataclass
class Conditions:
    """用户报告/否定的条件。空集合 = 条件未知（不是"没有"）。"""

    reported: frozenset[str] = frozenset()
    denied: frozenset[str] = frozenset()
    tags: frozenset[str] = frozenset()
    stage: str = ""
    notes: str = ""
    raw_reported: tuple[str, ...] = ()
    raw_denied: tuple[str, ...] = ()
    unknown_reported: tuple[str, ...] = ()

    def is_empty(self) -> bool:
        return not (self.reported or self.denied or self.tags or self.stage)

    def describe(self) -> str:
        if self.is_empty():
            return "条件未知（未报告任何试剂/步骤）"
        bits = []
        if self.reported:
            bits.append("已报告：" + "/".join(sorted(self.reported)))
        if self.denied:
            bits.append("明确否定：" + "/".join(sorted(self.denied)))
        if self.tags:
            bits.append("标签：" + "/".join(sorted(self.tags)))
        if self.stage:
            bits.append(f"阶段：{self.stage}")
        return "；".join(bits)


@dataclass
class Observed:
    """观测口径。kind=delta 时是"产物−起始物"的中性质量差。

    T04 扩充：
    - `mz_pair`：给两个实测峰（value、value2），由引擎求中性质量差并**传播误差**，
      归一后与 `delta` 完全同路径；`pair_derivation` 记录换算过程，报告里必显示。
    - `fixed_charge` 参与 m/z → 中性质量的换算：固定 +1 的季铵肽 [M]⁺ 不能按 [M+H]⁺ 反算。
    """

    kind: str
    value: float
    z: int = 0
    adduct: str = "H"
    tol_da: float | None = None
    tol_ppm: float | None = None
    value2: float | None = None
    note: str = ""
    # 归一化后的等价观测：统一换算成"中性质量差"或"中性质量"
    neutral: float | None = None
    # 双峰求差的换算说明（kind=mz 时为空）
    pair_derivation: str = ""
    # 双峰求差传播后的中性质量差误差（±Da）
    pair_error_da: float | None = None

    @property
    def effective_tol_da(self) -> float:
        if self.tol_da is not None:
            return float(self.tol_da)
        if self.pair_error_da is not None:
            return float(self.pair_error_da)
        if self.tol_ppm is not None:
            base = abs(self.value) if self.kind == "delta" else abs(self.value)
            return abs(base) * float(self.tol_ppm) * 1e-6
        # 默认：delta 用 0.02 Da（低分辨边界的一半），绝对质量用 10 ppm
        return 0.02 if self.kind == "delta" else 0.0

    def describe(self) -> str:
        if self.kind == "delta":
            s = f"ΔM = {self.value:+.6f} Da（容差 ±{self.effective_tol_da:.4f} Da）"
            if self.pair_derivation:
                s = f"ΔM = {self.value:+.6f} Da（{self.pair_derivation}；容差 ±{self.effective_tol_da:.4f} Da）"
        elif self.kind == "neutral_mass":
            ppm = self.tol_ppm if self.tol_ppm is not None else 10.0
            s = f"中性单同位素质量 M = {self.value:.6f} Da（容差 {ppm:g} ppm）"
        elif self.kind == "mz":
            s = f"m/z = {self.value:.6f}，z = {self.z}（加合 {self.adduct}，容差 {self.tol_ppm or 10:g} ppm）"
        else:  # pragma: no cover
            s = f"{self.kind}={self.value}"
        if self.note:
            s += f"；{self.note}"
        return s


@dataclass
class CaseInput:
    case_id: str
    state: PeptideState
    conditions: Conditions
    observed: Observed
    notes: str = ""
    fixed_charge: int = 0
    warnings: list[dict] = field(default_factory=list)


def parse_observed(raw: Any, fixed_charge: int = 0) -> Observed:
    if not isinstance(raw, dict):
        raise InputError("E_BAD_TYPE", "observed", "observed 必须是对象")
    _check_fields(raw, OBS_FIELDS, "observed")
    from .ions import ADDUCTS, neutral_from_ion
    from .formula import neutral_from_mz

    kind = str(raw.get("kind", "")).strip().lower()
    if kind in {"delta", "dm", "mass_shift", "质量差"}:
        kind = "delta"
    elif kind in {"neutral_mass", "mass", "m", "中性质量"}:
        kind = "neutral_mass"
    elif kind in {"mz", "m/z"}:
        kind = "mz"
    elif kind in {"mz_pair", "mz-pair", "双峰", "两个峰"}:
        kind = "mz_pair"
    else:
        raise InputError("E_BAD_OBSERVED_KIND", "observed.kind", f"未知观测类型 {raw.get('kind')!r}",
                         "用 'delta'、'neutral_mass'、'mz' 或 'mz_pair'（两个实测峰求差）")
    if "value" not in raw:
        raise InputError("E_MISSING_FIELD", "observed.value", "缺少观测值", "例如 {'kind':'delta','value':51.101428}")
    try:
        value = float(raw["value"])
    except (TypeError, ValueError) as exc:
        raise InputError("E_BAD_NUMBER", "observed.value", f"观测值不是数字：{raw['value']!r}") from exc

    z = int(raw.get("z", 0) or 0)
    z2 = int(raw.get("z2", 0) or 0)
    adduct = str(raw.get("adduct", "H") or "H")
    if adduct not in ADDUCTS:
        raise InputError("E_BAD_ADDUCT", "observed.adduct", f"未登记的加合离子 {adduct!r}",
                         f"可用：{', '.join(sorted(ADDUCTS))}（固定电荷肽用 fixed_charge 字段另说）")
    tol_da = raw.get("tol_da")
    tol_ppm = raw.get("tol_ppm")
    if tol_da is not None and float(tol_da) < 0:
        raise InputError("E_BAD_TOLERANCE", "observed.tol_da", "容差不能为负")
    if tol_ppm is not None and float(tol_ppm) < 0:
        raise InputError("E_BAD_TOLERANCE", "observed.tol_ppm", "容差不能为负")

    obs = Observed(
        kind=kind,
        value=value,
        z=z,
        adduct=adduct,
        tol_da=None if tol_da is None else float(tol_da),
        tol_ppm=None if tol_ppm is None else float(tol_ppm),
        value2=None if raw.get("value2") is None else float(raw["value2"]),
        note=str(raw.get("note", "") or ""),
    )

    if obs.kind == "mz":
        if obs.z == 0:
            raise InputError("E_MISSING_FIELD", "observed.z", "m/z 观测必须给电荷 z", "例如 {'kind':'mz','value':412.7,'z':2}")
        obs.neutral = neutral_from_ion(obs.value, obs.z, obs.adduct, fixed_charge)
    elif obs.kind == "neutral_mass":
        obs.neutral = obs.value
    elif obs.kind == "delta":
        obs.neutral = None  # delta 本身就是中性质量差
    else:  # mz_pair：两个实测峰求差，误差传播并核对电荷/加合
        if obs.value2 is None:
            raise InputError("E_MISSING_FIELD", "observed.value2", "mz_pair 需要第二个峰的 m/z", "value 与 value2 都填 m/z")
        if obs.z == 0:
            raise InputError("E_MISSING_FIELD", "observed.z", "mz_pair 需要给电荷 z（两峰相同）", "两个峰同一电荷态时只填 z")
        z2 = z2 or obs.z
        n1 = neutral_from_ion(obs.value, obs.z, obs.adduct, fixed_charge)
        n2 = neutral_from_ion(obs.value2, z2, obs.adduct, fixed_charge)
        obs.neutral = n2 - n1
        obs.value = n2 - n1
        obs.kind = "delta"
        # 误差传播：每个峰的 m/z 误差 → 中性质量误差 = m/z 误差 × z，方和根
        ppm = float(obs.tol_ppm) if obs.tol_ppm is not None else 10.0
        err1 = abs(raw["value"]) * ppm * 1e-6 * abs(obs.z)
        err2 = abs(raw["value2"]) * ppm * 1e-6 * abs(z2)
        obs.pair_error_da = (err1 ** 2 + err2 ** 2) ** 0.5
        obs.pair_derivation = (
            f"双峰求差：峰1 m/z {raw['value']:.6f} @ z={obs.z} → 中性 {n1:.6f}；"
            f"峰2 m/z {raw['value2']:.6f} @ z={z2} → 中性 {n2:.6f}；"
            f"差值 {n2 - n1:+.6f} Da；误差按 {ppm:g} ppm 传播为 ±{obs.pair_error_da:.6f} Da"
        )
        if obs.tol_da is None:
            obs.tol_ppm = None   # 已在 pair_error_da 里体现，避免重复放大
    return obs


def parse_conditions(raw: Any) -> Conditions:
    if raw is None:
        return Conditions()
    if not isinstance(raw, dict):
        raise InputError("E_BAD_TYPE", "conditions", "conditions 必须是对象")
    _check_fields(raw, COND_FIELDS, "conditions")

    reported_raw: list[str] = []
    denied_raw: list[str] = []
    for key, bucket in (("reported_reagents", reported_raw), ("denied_reagents", denied_raw)):
        vals = raw.get(key, []) or []
        if isinstance(vals, str):
            vals = [vals]
        if not isinstance(vals, list):
            raise InputError("E_BAD_TYPE", f"conditions.{key}", "应为字符串列表")
        bucket.extend(str(v) for v in vals)

    reported: set[str] = set()
    unknown_reported: list[str] = []
    for blob in reported_raw:
        for piece in split_reagent_blob(blob):
            canon = canon_reagent(piece) or _keep_unknown_reagent(piece)
            if canon:
                reported.add(canon)
            if piece.strip().lower() not in TAG_IS_MODIFICATION and canon_reagent(piece) is None:
                unknown_reported.append(piece)
    denied: set[str] = set()
    for blob in denied_raw:
        for piece in split_reagent_blob(blob):
            canon = canon_reagent(piece) or _keep_unknown_reagent(piece)
            if canon:
                denied.add(canon)
            if canon_reagent(piece) is None:
                unknown_reported.append(piece)

    tags_raw = raw.get("tags", []) or []
    if isinstance(tags_raw, str):
        tags_raw = [tags_raw]
    tags: set[str] = set()
    for t in tags_raw:
        key = str(t).strip().lower()
        canon = canon_reagent(key)
        if canon:
            # 标签里写试剂名也算报告了该试剂
            reported.add(canon)
            tags.add(canon.lower())
        else:
            tags.add(key)

    stage = str(raw.get("stage", "") or "")
    notes = str(raw.get("notes", "") or "")
    return Conditions(
        reported=frozenset(reported),
        denied=frozenset(denied),
        tags=frozenset(tags),
        stage=stage,
        notes=notes,
        raw_reported=tuple(reported_raw),
        raw_denied=tuple(denied_raw),
        unknown_reported=tuple(unknown_reported),
    )


def parse_case(raw: Any, lib: ComponentLibrary) -> CaseInput:
    """解析整条病例。严格白名单 + 明确报错。"""
    if not isinstance(raw, dict):
        raise InputError("E_BAD_TYPE", "case", "顶层应为对象（dict）")
    _check_fields(raw, TOP_FIELDS, "case")
    if "sequence" not in raw:
        raise InputError("E_MISSING_FIELD", "sequence", "缺少序列", "用一字母序列，如 'GC'、'Ac-DYKDDDDK-NH2'")

    mods = raw.get("modifications", []) or []
    if isinstance(mods, dict):
        raise InputError("E_BAD_TYPE", "modifications", "应为列表", '形如 [{"site":"N","group":"Fmoc"}]')
    if not isinstance(mods, list):
        raise InputError("E_BAD_TYPE", "modifications", "应为列表")
    for i, m in enumerate(mods):
        if not isinstance(m, dict):
            raise InputError("E_BAD_TYPE", f"modifications[{i}]", "每个修饰应为对象")
        _check_fields(m, MOD_FIELDS, f"modifications[{i}]")

    warnings: list[dict] = []
    try:
        state = build_state(
            sequence=raw["sequence"],
            nterm=raw.get("nterm") or "free",
            cterm=raw.get("cterm") or "acid",
            modifications=mods,
            lib=lib,
        )
    except StateError as exc:
        raise InputError(exc.code, "structure", exc.message, exc.hint) from exc

    # 端基通过 nterm 字段给了，又用 modifications 给了 N 端基团 → 重复，直接报错
    if state.nterm_group and any(m.site == "N" for m in state.modifications):
        raise InputError(
            "E_DUPLICATE_NTERM",
            "nterm",
            f"N 端基团重复：nterm={state.nterm_group}，同时又有 N 端修饰 "
            + ",".join(m.group for m in state.modifications if m.site == "N"),
            "二者只写一处；推荐统一用 nterm 或统一用 modifications",
        )

    conditions = parse_conditions(raw.get("conditions"))
    if "observed" not in raw:
        raise InputError("E_MISSING_FIELD", "observed", "缺少观测（质量差或实测质量）", "例如 {'kind':'delta','value':51.101428}")

    fixed_charge = raw.get("fixed_charge", 0)
    if fixed_charge is None:
        fixed_charge = 0
    try:
        fixed_charge = int(fixed_charge)
    except (TypeError, ValueError) as exc:
        raise InputError("E_BAD_FIXED_CHARGE", "fixed_charge", f"固定电荷必须是整数，收到 {raw.get('fixed_charge')!r}") from exc
    if not -4 <= fixed_charge <= 4:
        raise InputError("E_BAD_FIXED_CHARGE", "fixed_charge", f"固定电荷 {fixed_charge} 超出合理范围（−4…+4）",
                         "只有永久带电基团（季铵、锍、吡啶鎓、TMPP 等）才写这个字段；质子化不算固定电荷")
    if fixed_charge:
        warnings.append({
            "code": "W_FIXED_CHARGE",
            "message": f"按永久固定电荷 {fixed_charge:+d} 处理：m/z ↔ 中性质量换算不再加质子（与 [M+H]+ 的差约为 1.0073/z Da）",
        })

    observed = parse_observed(raw.get("observed"), fixed_charge=fixed_charge)

    for piece in conditions.unknown_reported:
        warnings.append({
            "code": "W_UNKNOWN_REAGENT",
            "message": f"未识别的试剂/条件写法 {piece!r}，未参与条件判断（不会被当成'没有'）",
        })
    if conditions.denied and not conditions.reported:
        warnings.append({
            "code": "W_ONLY_NEGATIVES",
            "message": "只给了否定条件；据此可以排除依赖该试剂的规则，但无法确认任何规则的阳性条件",
        })
    if state.folded:
        warnings.append({"code": "W_FOLDED", "message": "；".join(state.folded)})

    return CaseInput(
        case_id=str(raw.get("case_id", "") or "case"),
        state=state,
        conditions=conditions,
        observed=observed,
        notes=str(raw.get("notes", "") or ""),
        fixed_charge=fixed_charge,
        warnings=warnings,
    )