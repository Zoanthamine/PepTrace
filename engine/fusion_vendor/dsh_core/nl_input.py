"""T08 自然语言入口（受限、可解释）。

目标：把中文实验记录里的几句话转成病例 dict，再交给 parse_case 走与 JSON 输入
完全相同的校验和推理路径——自然语言层不做任何化学判断，也不补全没说的条件。

支持的句式（其余一律报错并给出可用的写法，绝不猜）：
- 序列：连续的氨基酸单字母大写串（≥3 个），如 ``序列 MELK`` / ``MELK 肽``
- 质量差：``ΔM=+15.99 Da`` / ``多出 15.99 Da`` / ``少了 18.01 Da`` / ``质量增加 79.97 amu``
- 绝对质量：``实测 535.2676 Da`` / ``测得 519.3 Da``
- 前体离子：``m/z 520.2799``（可带 ``z=1`` 或 ``2+``）
- 试剂/条件：命中 REAGENT_ALIASES 的英文别名或中文名（哌啶、三氟乙酸、乙二硫醇…）
- 端基：``C 端酰胺`` / ``C端为酰胺``

明确不做：序列歧义消解、单位换算猜测、隐含条件补全、多候选质量解析。
"""

from __future__ import annotations

import re

from .schema import REAGENT_ALIASES, InputError

__all__ = ["NLParseError", "parse_text"]

_SEQ_RE = re.compile(r"(?<![A-Za-z])([ACDEFGHIKLMNPQRSTVWY]{3,})(?![A-Za-z])")
_DELTA_WORD_RE = re.compile(
    r"(ΔM|Δm|delta|质量差|多出|多了|增加|升高|少了|减少|降低|相差|差值)"
    r"[^0-9]{0,8}([0-9]+(?:\.[0-9]+)?)\s*(?:Da|da|amu|道尔顿)?",
    re.IGNORECASE,
)
_DELTA_SIGN_RE = re.compile(
    r"([+＋\-−])\s*([0-9]+(?:\.[0-9]+)?)\s*(?:Da|da|amu|道尔顿)",
    re.IGNORECASE,
)
_NEGATIVE_WORDS = {"少了", "减少", "降低"}
_REAGENT_UPPER = {k.upper() for k in REAGENT_ALIASES if k.isascii() and len(k) >= 2}
_ABS_RE = re.compile(r"(?:实测|测得|检测到|测定为|绝对质量)\s*[约为:：]?\s*([0-9]{2,5}\.[0-9]{2,6})\s*(Da|da)?")
_MZ_RE = re.compile(r"m/z\s*[=为:：]?\s*([0-9]+\.[0-9]+)", re.IGNORECASE)
_Z_RE = re.compile(r"(?:z\s*=\s*|\b)([1-4])\s*\+")
_CTERM_AMIDE_RE = re.compile(r"C\s*端\s*(?:为|是)?\s*酰胺")

_CN_REAGENTS = {
    "哌啶": "Piperidine",
    "三氟乙酸": "TFA",
    "TFA": "TFA",
    "乙二硫醇": "EDT",
    "EDT": "EDT",
    "三异丙基硅烷": "TIS",
    "TIS": "TIS",
    "二异丙基乙胺": "DIEA",
    "DIEA": "DIEA",
    "DMSO": "DMSO",
    "二甲基亚砜": "DMSO",
    "甲醇": "Methanol",
    "乙腈": "ACN",
    "ACN": "ACN",
    "二氯甲烷": "DCM",
    "DMF": "DMF",
    "NMP": "NMP",
    "MTBE": "MTBE",
    "肼": "Hydrazine",
    "胍": "Guanidine",
    "甲醛": "Formaldehyde",
    "丙酮": "Acetone",
    "苯甲硫醚": "Thioanisole",
    "苯酚": "Phenol",
    "氢氟酸": "HF",
    "HF": "HF",
    "甲酸": "FormicAcid",
    "乙酸": "AceticAcid",
}


class NLParseError(InputError):
    """自然语言输入无法可靠解析。错误码固定为 E_NL_PARSE。"""

    def __init__(self, message: str, hint: str):
        super().__init__("E_NL_PARSE", "text", message, hint)


def _find_reagents(text: str) -> list[str]:
    found: list[str] = []
    upper = text.upper()
    for alias, canonical in _CN_REAGENTS.items():
        token = alias if alias in text else (alias if alias.upper() in upper else None)
        if token and canonical not in found:
            found.append(canonical)
    return found


def _sentences_split(text: str) -> list[str]:
    return [s for s in re.split(r"[的，,。；;、\s]+", text) if s]


def parse_text(text: str, *, case_id: str | None = None) -> dict:
    """把一句话（或几行记录）解析成病例 dict。失败抛 NLParseError，不返回半成品。"""
    if not text or not text.strip():
        raise NLParseError("输入为空", "直接写「序列 MELK，ΔM=+15.99 Da」这样的句子")

    seq_matches = [
        m.group(1) for m in _SEQ_RE.finditer(text)
        if m.group(1).upper() not in _REAGENT_UPPER
    ]
    sequences = list(dict.fromkeys(seq_matches))
    if not sequences:
        raise NLParseError(
            "没有识别到氨基酸序列",
            "用单字母大写写序列，如「序列 MELK」；不要用三字母或中文残基名",
        )
    if len(sequences) > 1:
        raise NLParseError(
            f"识别到多个候选序列：{sequences}",
            "一句话里只留一个序列；若要比较多个肽请分别建病例",
        )
    sequence = sequences[0]

    observed: dict[str, object] = {}
    notes: list[str] = []

    mz = _MZ_RE.search(text)
    if mz:
        z = 1
        zm = _Z_RE.search(text)
        if zm:
            z = int(zm.group(1))
        observed = {"kind": "mz", "value": float(mz.group(1)), "z": z, "adduct": "H"}
        obs_vals = [m.group(1) for m in re.finditer(r"([0-9]+\.[0-9]+)", text)]
        others = [v for v in obs_vals if v != mz.group(1)]
        if others:
            notes.append(f"文本里另有数值 {others}，未当作第二个峰；若是双峰求差请用 JSON 的 mz_pair")
        notes.append("前体按 [M+H]+ 解释；若实际是钠/钾加合请在 JSON 里写 adduct")
    else:
        ab = _ABS_RE.search(text)
        word = _DELTA_WORD_RE.search(text)
        sign_m = _DELTA_SIGN_RE.search(text)
        if word and not (ab and ab.start() < word.start()):
            key = (word.group(1) or "").lower()
            sign = -1.0 if key in _NEGATIVE_WORDS else 1.0
            before = text[max(0, word.start(2) - 3):word.start(2)]
            if "-" in before or "−" in before:
                sign = -1.0
            elif "+" in before or "＋" in before:
                sign = 1.0
            observed = {"kind": "delta", "value": sign * float(word.group(2))}
            if key not in {"δm", "δm", "delta", "质量差", "δm"} and word.group(1) not in {"ΔM", "Δm"}:
                notes.append(f"质量变化方向取自「{word.group(1)}」；建议直接写 +/− 号")
        elif sign_m:
            sign = -1.0 if sign_m.group(1) in {"-", "−"} else 1.0
            observed = {"kind": "delta", "value": sign * float(sign_m.group(2))}
        elif ab:
            observed = {"kind": "neutral_mass", "value": float(ab.group(1))}
            notes.append("按绝对质量解释，结果为实测−目标")
        else:
            raise NLParseError(
                "没有识别到质量信息",
                "写清楚是质量差（ΔM=+15.99 Da）还是绝对质量（实测 535.2676 Da）或前体 m/z；"
                "质量差请带 +/− 号或「多出/少了」这类量词",
            )

    denials = re.findall(r"(?:没有|未|不使用|不加|不含)(?:使用)?\s*([\u4e00-\u9fa5A-Za-z]{2,8})", text)
    denied: list[str] = []
    for d in denials:
        denied.extend(_find_reagents(d))
    denied = list(dict.fromkeys(denied))
    reagents = [r for r in _find_reagents(text) if r not in denied]
    conditions: dict[str, object] = {}
    if reagents:
        conditions["reported_reagents"] = reagents
    if denied:
        conditions["denied_reagents"] = denied
        notes.append(f"识别到阴性描述（{denied}），已作为反证字段，不静默忽略")

    raw: dict[str, object] = {"sequence": sequence, "observed": observed}
    if case_id:
        raw["case_id"] = case_id
    if conditions:
        raw["conditions"] = conditions
    if _CTERM_AMIDE_RE.search(text):
        raw["cterm"] = "amide"
        notes.append("识别到 C 端酰胺，已按 cterm=amide 归一到肽质量")
    if notes:
        raw["notes"] = "；".join(notes)
    return raw