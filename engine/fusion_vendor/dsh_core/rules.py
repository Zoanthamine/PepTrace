"""规则库加载与规范化。

只读两份既有资产（不改写）：
- `规则库_总.json`：382 条主规则（来自全书 14 章）
- `SUPPLEMENTS.json`：63 条编辑补充条目（origin=supplement，不得当原书数据引用）

加载时做两件"自证"工作：
1. 用本包独立原子量表重算 `delta_formula`，与 `delta_mass_mono_da` 比对；
   不一致就记入 `problems`（**不静默取其一**，报告里单列"数据自检"）。
2. 把中文 `nominal_mass_shift_book` 解析成整数与"倍数"标记，供低分辨检索。
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from .components import ComponentLibrary
from .formula import Formula, formula_str

_NOMINAL_RE = re.compile(r"([+\-−–])\s*([0-9]+(?:\.[0-9]+)?)")


@dataclass
class Rule:
    rule_id: str
    name_cn: str
    name_en: str
    chapter: int
    section: str
    pdf_pages: tuple[int, ...]
    print_pages: tuple[int, ...]
    figures: tuple[str, ...]
    reference_numbers: tuple[str, ...]
    origin: str  # book | supplement
    stage: tuple[str, ...]
    conditions: str
    products: str
    susceptible: tuple[str, ...]
    mechanism: str
    mass_status: str
    mass_expression: str
    mass_reference: str
    diagnostic_evidence: str
    competing: str
    mitigation: str
    certainty: str
    # 规范化后的质量信息
    delta: Formula | None = None
    delta_mass: float | None = None
    nominal: float | None = None
    nominal_multiples: bool = False
    # 书里的原始整数口径文字（如 '−42 amu'、'+80 amu 及其倍数'），必须原样保留：
    # 低分辨检索与报告都要能引用书中的整数说法，不能被规范化后丢掉。
    nominal_book: str = ""
    polymer_relation: str = ""
    # 条件与残基要求（来自 conditions_map.json）
    required_any: tuple[str, ...] = ()
    contextual: tuple[str, ...] = ()
    named_residues: tuple[str, ...] = ()

    @property
    def has_exact_mass(self) -> bool:
        return self.delta_mass is not None

    @property
    def is_low_res_only(self) -> bool:
        return self.delta_mass is None and self.nominal is not None

    def mass_label(self) -> str:
        if self.delta_mass is not None:
            return f"{self.delta_mass:+.6f} Da"
        if self.nominal is not None:
            return f"{self.nominal:+.0f} amu（仅整数）"
        return "无单一质量差"


@dataclass
class RuleSet:
    rules: dict[str, Rule]
    problems: list[dict] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.rules)

    def by_chapter(self) -> dict[int, list[Rule]]:
        out: dict[int, list[Rule]] = {}
        for r in self.rules.values():
            out.setdefault(r.chapter, []).append(r)
        return out

    def mass_matchable(self) -> list[Rule]:
        """能参与数值比对的规则：有精确质量差，或有整数 amu。"""
        return [r for r in self.rules.values() if r.has_exact_mass or r.is_low_res_only]

    def stats(self) -> dict:
        from collections import Counter

        return {
            "total": len(self.rules),
            "book": sum(1 for r in self.rules.values() if r.origin == "book"),
            "supplement": sum(1 for r in self.rules.values() if r.origin == "supplement"),
            "exact_mass": sum(1 for r in self.rules.values() if r.has_exact_mass),
            "nominal_only": sum(1 for r in self.rules.values() if r.is_low_res_only),
            "no_mass": sum(1 for r in self.rules.values() if not (r.has_exact_mass or r.is_low_res_only)),
            "mass_status": dict(Counter(r.mass_status for r in self.rules.values())),
            "data_problems": len(self.problems),
        }


# 多聚/二聚关系：[2M−18] 这类写的是「产物总质量 = 2×单体 − 18」，不是单体的质量差。
# 把它当成 −18 Da 的单事件会把整库检索搞错（T06 独立穷举对照抓到的真实缺陷）。
_POLYMER_RE = re.compile(r"\[[2-9]\s*M|\[\s*[2-9]\s*×\s*M|多聚体|二聚体")


def polymer_relation(text) -> str:
    """若文本描述的是多聚/二聚平衡关系，返回该描述原文；否则返回空串。"""
    if not text:
        return ""
    s = str(text)
    m = _POLYMER_RE.search(s)
    return s if m else ""


def _parse_nominal(text) -> tuple[float | None, bool]:
    if not text:
        return None, False
    s = str(text)
    if polymer_relation(s):
        # [2M−18] 是「产物总质量」关系，不能降级成单体净质量差；交给关系层单独说明
        return None, False
    m = _NOMINAL_RE.search(s)
    if not m:
        return None, False
    sign = -1.0 if m.group(1) in "-−–" else 1.0
    val = sign * float(m.group(2))
    multiples = ("倍数" in s) or ("×" in s) or ("n×" in s)
    return val, multiples


def _txt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return "；".join(str(x) for x in v)
    return str(v)


def _tup(v) -> tuple[str, ...]:
    if v is None:
        return ()
    if isinstance(v, (list, tuple)):
        return tuple(str(x) for x in v)
    return (str(v),)


def _ints(v) -> tuple[int, ...]:
    out = []
    for x in _tup(v):
        try:
            out.append(int(float(x)))
        except (TypeError, ValueError):
            continue
    return tuple(out)


def load_rules(lib: ComponentLibrary, kb_root: Path | None = None) -> RuleSet:
    root = kb_root or Path(__file__).resolve().parents[2]
    book_path = root / "规则库_总.json"
    sup_path = root / "SUPPLEMENTS.json"
    cond_path = Path(__file__).resolve().parent / "conditions_map.json"
    if not book_path.exists():
        raise FileNotFoundError(f"找不到规则库：{book_path}")

    raw_rules = json.loads(book_path.read_text(encoding="utf-8"))["rules"]
    raw_sups = json.loads(sup_path.read_text(encoding="utf-8"))["items"] if sup_path.exists() else []
    cond = json.loads(cond_path.read_text(encoding="utf-8"))["entries"] if cond_path.exists() else {}

    rules: dict[str, Rule] = {}
    problems: list[dict] = []
    for raw in list(raw_rules) + list(raw_sups):
        rid = raw["rule_id"]
        src = raw.get("source") or {}
        origin = raw.get("origin", "book")
        mask = raw.get("mass_status", "")
        df = raw.get("delta_formula")
        delta: Formula | None = None
        if isinstance(df, dict) and df:
            try:
                delta = Formula.of(df)
            except Exception as exc:  # noqa: BLE001
                problems.append({"rule_id": rid, "kind": "bad_delta_formula", "detail": f"{df} → {exc}"})
        stored = raw.get("delta_mass_mono_da")
        stored = float(stored) if stored is not None else None
        if delta is not None and stored is not None:
            mine = delta.mass()
            if abs(mine - stored) > 1e-3:
                problems.append({
                    "rule_id": rid,
                    "kind": "delta_mass_mismatch",
                    "detail": f"delta_formula '{delta.name()}' 重算 {mine:.6f}，库中记录 {stored:.6f}，差 {mine - stored:+.6f}",
                })
        nominal_book = _txt(raw.get("nominal_mass_shift_book"))
        nominal, multiples = _parse_nominal(nominal_book)
        c = cond.get(rid, {})
        rules[rid] = Rule(
            rule_id=rid,
            name_cn=_txt(raw.get("name_cn")),
            name_en=_txt(raw.get("name_en")),
            chapter=int(src.get("chapter") or raw.get("_chapter") or 0),
            section=_txt(src.get("section")),
            pdf_pages=_ints(src.get("pdf_pages")),
            print_pages=_ints(src.get("print_pages")),
            figures=_tup(src.get("figures")),
            reference_numbers=_tup(src.get("reference_numbers")),
            origin=origin,
            stage=_tup(raw.get("stage")),
            conditions=_txt(raw.get("conditions")),
            products=_txt(raw.get("products")),
            susceptible=_tup(raw.get("susceptible_sequences")),
            mechanism=_txt(raw.get("mechanism_summary")),
            mass_status=mask,
            mass_expression=_txt(raw.get("mass_expression")),
            mass_reference=_txt(raw.get("mass_reference")),
            diagnostic_evidence=_txt(raw.get("diagnostic_evidence")),
            competing=_txt(raw.get("competing_explanations")),
            mitigation=_txt(raw.get("mitigation")),
            certainty=_txt(raw.get("certainty_notes")),
            delta=delta,
            delta_mass=stored if stored is not None else (delta.mass() if delta else None),
            nominal=nominal,
            nominal_multiples=multiples,
            nominal_book=nominal_book,
            polymer_relation=polymer_relation(nominal_book) or polymer_relation(raw.get("mass_expression")),
            required_any=tuple(c.get("required_any") or ()),
            contextual=tuple(c.get("contextual") or ()),
            named_residues=tuple(c.get("named_residues") or ()),
        )
    return RuleSet(rules=rules, problems=problems)


def rule_brief(r: Rule) -> dict:
    """报告用精简结构（不含长文本）。"""
    return {
        "rule_id": r.rule_id,
        "name_cn": r.name_cn,
        "chapter": r.chapter,
        "section": r.section,
        "print_pages": list(r.print_pages),
        "origin": r.origin,
        "mass_status": r.mass_status,
        "mass_label": r.mass_label(),
        "delta_formula": r.delta.name() if r.delta else None,
    }


def formula_mass_text(f: Formula | None) -> str:
    return formula_str(f.as_dict()) if f else ""