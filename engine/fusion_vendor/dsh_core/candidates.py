"""候选生成与证据判定：本项目的核心。

判定模型（一条候选的四件事，缺一不可）
--------------------------------------
1. **质量**：数值上是否吻合（exact / nominal）。不吻合 = 不进候选。
2. **适用性**：序列是否具备该规则要求的残基。缺残基 = 该候选不适用（排除，但报告里写明）。
3. **条件**：用户报告/否定的条件与规则要求的关系 →
   satisfied > weak（仅旁证）> unknown（**保留**）> violated（**排除**）。
4. **位点**：候选位点（残基位置），只说"可能在哪"，不假装能定结构。

关键政策（对应执行计划第 4 节六个失败案例）
- 条件未知**不删**候选，只降排序：`unknown` 仍然列出，并明确标注"条件待确认"。
- 阴性证据**参与**判断：用户明确说没用某试剂，依赖该试剂的规则转入"已排除"，并给出理由。
- 条件吻合但质量不符的规则**不静默丢弃**，作为反证单列，帮助用户判断"这个质量不像书上哪条"。
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .components import ComponentLibrary
from .formula import Formula
from .ions import mz_of
from .isotopes import envelope
from .matching import MassMatch, group_same_delta, match_absolute, match_delta
from .rules import Rule, RuleSet, rule_brief
from .schema import CaseInput, canon_reagent
from .state import PeptideState

COND_RANK = {"satisfied": 3, "weak": 2, "unknown": 1, "violated": 0}
COND_LABEL = {
    "satisfied": "条件支持",
    "weak": "条件旁证",
    "unknown": "条件未知（保留）",
    "violated": "条件冲突（排除）",
}


@dataclass
class Candidate:
    rule: Rule
    match: MassMatch
    applicability: str          # applicable | residue_missing
    applicability_reason: str
    condition_state: str
    condition_reason: str
    sites: list[dict]
    verdict: str                # retained | excluded
    excluded_reason: str = ""
    expected_product_mass: float | None = None
    # T04：产物元素组成（规则给了净式才有）、理论 m/z（m/z 观测才有）、同位素包络
    product_formula: Formula | None = None
    expected_mz: float | None = None
    isotope_lines: tuple[dict, ...] = ()

    @property
    def rule_id(self) -> str:
        return self.rule.rule_id

    def rank_key(self) -> tuple:
        grade = 0 if self.match.grade.startswith("exact") else 1
        err = abs(self.match.error_ppm) if self.match.error_ppm is not None else abs(self.match.error_da) * 1e3
        return (grade, -COND_RANK[self.condition_state], err)

    def as_row(self) -> dict:
        return {
            "rule_id": self.rule_id,
            "name_cn": self.rule.name_cn,
            "chapter": self.rule.chapter,
            "print_pages": list(self.rule.print_pages),
            "origin": self.rule.origin,
            "predicted_delta_da": round(self.match.predicted, 6),
            "error_da": round(self.match.error_da, 6),
            "error_ppm": None if self.match.error_ppm is None else round(self.match.error_ppm, 2),
            "mass_grade": self.match.grade,
            "match_note": self._match_note(),
            "mass_status": self.rule.mass_status,
            "condition_state": self.condition_state,
            "condition_label": COND_LABEL[self.condition_state],
            "condition_reason": self.condition_reason,
            "applicability": self.applicability,
            "sites": self.sites,
            "expected_product_mass": None if self.expected_product_mass is None else round(self.expected_product_mass, 6),
            "product_formula": None if self.product_formula is None else self.product_formula.name(),
            "expected_mz": None if self.expected_mz is None else round(self.expected_mz, 6),
            "isotope_lines": list(self.isotope_lines),
            "verdict": self.verdict,
            "excluded_reason": self.excluded_reason,
        }

    def _match_note(self) -> str:
        if self.match.grade == "nominal":
            return "按整数 amu 匹配（该规则只有整数质量，不能当精确匹配）"
        if self.match.grade == "exact_multiple":
            return f"按 {self.match.unit_ratio:g} 个重复单元的精确质量匹配（该规则标注「及其倍数」）"
        if self.match.grade == "nominal_multiple":
            return f"按 {self.match.unit_ratio:g} 个重复单元匹配（该规则标注'及其倍数'）"
        return "精确质量差匹配"


@dataclass
class Diagnosis:
    case: CaseInput
    target_formula_name: str
    target_mass: float
    observed_desc: str
    observed_delta: float | None
    tol_da: float
    candidates: list[Candidate] = field(default_factory=list)
    excluded: list[Candidate] = field(default_factory=list)
    counter_evidence: list[dict] = field(default_factory=list)
    groups: list[dict] = field(default_factory=list)
    data_problems: list[dict] = field(default_factory=list)
    stats: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def retained_ids(self) -> list[str]:
        return [c.rule_id for c in self.candidates]

    def by_id(self, rid: str) -> Candidate | None:
        for c in self.candidates + self.excluded:
            if c.rule_id == rid:
                return c
        return None


def _canon_set(values) -> set[str]:
    """把规则派生条件与用户输入都归一到同一套规范名。

    踩过的真 bug：CH07-R030 的 required_any 是从正文抽出的 'formaldehyde'（小写），
    而用户否定 'Formaldehyde' 被归一成规范名，两边直接求交会漏掉阴性证据——
    即"用户说没用甲醛，规则却仍被当作未知条件保留"。
    """
    out: set[str] = set()
    for v in values or ():
        c = canon_reagent(str(v)) or str(v)
        out.add(c)
    return out


def _eval_conditions(rule: Rule, case: CaseInput) -> tuple[str, str]:
    cond = case.conditions
    req = _canon_set(rule.required_any)
    reported = _canon_set(cond.reported)
    denied_set = _canon_set(cond.denied)
    contextual = _canon_set(rule.contextual)
    denied = req & denied_set
    if denied:
        return "violated", (
            f"你明确否定使用 {'/'.join(sorted(denied))}，而该规则要求 {'/'.join(sorted(req))} 之一；"
            "按阴性证据将该候选排除（不是忽略你的说明）"
        )
    hit = req & reported
    if hit:
        rest = req - hit
        extra = f"（该规则还提到 {'/'.join(sorted(rest))}，你未报告；有其一即可支持）" if rest else ""
        return "satisfied", f"你报告了 {'/'.join(sorted(hit))}，与该规则条件一致{extra}"
    if req:
        return "unknown", (
            f"该规则需要 {'/'.join(sorted(req))} 之一，你既未报告也未否定；**保留候选、条件待确认**"
        )
    chit = contextual & (reported | denied_set)
    if chit - denied_set:
        return "weak", f"该规则正文提到 {'/'.join(sorted(chit - denied_set))}，与你的报告一致；无权威条件要求，仅作旁证"
    if cond.is_empty():
        return "unknown", "未提供任何条件；按'条件未知保留'政策保留，不因缺信息删候选"
    return "unknown", "你报告的条件与该规则无交集，既不支持也不否定；保留候选"


def _eval_residues(rule: Rule, state: PeptideState) -> tuple[str, str, list[dict]]:
    seq = state.sequence
    if not rule.named_residues:
        return "applicable", "规则未限定特定残基", []
    present = [aa for aa in rule.named_residues if aa.upper() in seq.upper()]
    if not present:
        return (
            "residue_missing",
            f"规则涉及 {'/'.join(sorted(rule.named_residues))}，而序列 {'/'.join(sorted(set(seq)))} 一个都不含",
            [],
        )
    sites: list[dict] = []
    if len(present) <= 4:
        for aa in present:
            if aa == "pGlu":
                continue
            pos = state.positions_of(aa)
            if pos:
                sites.append({"residue": aa, "positions": pos})
    return "applicable", f"序列含该规则涉及的残基 {'/'.join(sorted(present))}", sites


def diagnose(case: CaseInput, ruleset: RuleSet, lib: ComponentLibrary) -> Diagnosis:
    state = case.state
    obs = case.observed
    target_formula = state.formula(lib)
    target_mass = target_formula.mass()

    if obs.kind == "delta":
        observed_delta = obs.value
        tol_da = obs.effective_tol_da
        matches = match_delta(ruleset.mass_matchable(), observed_delta, tol_da)
    else:
        assert obs.neutral is not None, "绝对质量观测必须已归一为中性质量"
        ppm = obs.tol_ppm if obs.tol_ppm is not None else 10.0
        tol_da = abs(obs.neutral) * float(ppm) * 1e-6
        observed_delta = obs.neutral - target_mass
        matches = match_absolute(ruleset.mass_matchable(), target_mass, obs.neutral, tol_da)

    diag = Diagnosis(
        case=case,
        target_formula_name=target_formula.name(),
        target_mass=target_mass,
        observed_desc=obs.describe(),
        observed_delta=observed_delta,
        tol_da=tol_da,
        data_problems=ruleset.problems,
    )

    retained: list[Candidate] = []
    excluded: list[Candidate] = []
    for m in matches:
        rule = m.rule
        app, app_reason, sites = _eval_residues(rule, state)
        cstate, creason = _eval_conditions(rule, case)
        verdict = "retained"
        excluded_reason = ""
        if app == "residue_missing":
            verdict, excluded_reason = "excluded", f"序列不具备该规则的残基前提：{app_reason}"
        elif cstate == "violated":
            verdict, excluded_reason = "excluded", creason
        # T04：产物组成 / 理论 m/z / 同位素包络
        product_formula = None
        if rule.delta is not None:
            product_formula = target_formula + rule.delta * int(m.unit_ratio or 1)
        expected_mz = None
        isotope_lines: tuple[dict, ...] = ()
        if obs.kind == "mz" and obs.z:
            product_mass = (product_formula.mass() if product_formula is not None
                            else target_mass + m.predicted)
            expected_mz = mz_of(product_mass, obs.z, obs.adduct, fixed_charge=case.fixed_charge)
            if product_formula is not None:
                env = envelope(product_formula, max_offset=2)
                isotope_lines = tuple({
                    "offset": row["offset"],
                    "mz": round(expected_mz + row["mass_shift"] / abs(obs.z), 6),
                    "rel_intensity": row["rel_intensity"],
                } for row in env)
        cand = Candidate(
            rule=rule,
            match=m,
            applicability=app,
            applicability_reason=app_reason,
            condition_state=cstate,
            condition_reason=creason,
            sites=sites,
            verdict=verdict,
            excluded_reason=excluded_reason,
            expected_product_mass=target_mass + m.predicted,
            product_formula=product_formula,
            expected_mz=expected_mz,
            isotope_lines=isotope_lines,
        )
        (retained if verdict == "retained" else excluded).append(cand)

    retained.sort(key=lambda c: c.rank_key())
    excluded.sort(key=lambda c: (c.condition_state != "violated", abs(c.match.error_da)))
    diag.candidates = retained
    diag.excluded = excluded
    diag.groups = group_same_delta([c.match for c in retained])

    # ---- 反证：条件吻合、序列适用，但质量对不上 ----
    matched_ids = {c.rule_id for c in retained + excluded}
    contra: list[dict] = []
    if observed_delta is not None:
        for rule in ruleset.rules.values():
            if rule.rule_id in matched_ids:
                continue
            if rule.delta_mass is None and rule.nominal is None:
                continue
            app, _, _ = _eval_residues(rule, state)
            if app != "applicable":
                continue
            cstate, _ = _eval_conditions(rule, case)
            if cstate not in {"satisfied", "weak"}:
                continue
            base = rule.delta_mass if rule.delta_mass is not None else rule.nominal
            assert base is not None
            contra.append({
                "rule_id": rule.rule_id,
                "name_cn": rule.name_cn,
                "chapter": rule.chapter,
                "predicted_delta_da": round(base, 6),
                "gap_da": round(observed_delta - base, 6),
                "why": f"条件与序列都吻合（{COND_LABEL[cstate]}），但观测差 {observed_delta:+.6f} 与预测 {base:+.6f} 相差 {observed_delta - base:+.6f} Da",
            })
    contra.sort(key=lambda x: abs(x["gap_da"]))
    diag.counter_evidence = contra[:8]

    diag.stats = {
        "rule_total": len(ruleset),
        "mass_matchable": len(ruleset.mass_matchable()),
        "mass_matched": len(matches),
        "retained": len(retained),
        "excluded": len(excluded),
        "excluded_condition_conflict": sum(1 for c in excluded if c.condition_state == "violated"),
        "excluded_residue_missing": sum(1 for c in excluded if c.applicability == "residue_missing"),
        "condition_unknown_retained": sum(1 for c in retained if c.condition_state == "unknown"),
        "condition_satisfied_retained": sum(1 for c in retained if c.condition_state == "satisfied"),
        "same_delta_groups": len(diag.groups),
        "counter_evidence": len(diag.counter_evidence),
        "observed_delta": None if observed_delta is None else round(observed_delta, 6),
        "tolerance_da": round(tol_da, 6),
    }
    if not retained:
        diag.notes.append(
            "在声明容差内没有质量吻合的规则。可能原因：序列/端基/修饰输入与真实起始物不一致；"
            "电荷或加合离子判断有误；观测峰不是产物峰；或该副反应不在本书范围内。"
        )
    if diag.stats["condition_unknown_retained"]:
        diag.notes.append(
            f"有 {diag.stats['condition_unknown_retained']} 条候选因**条件未知**而保留。"
            "补上试剂/步骤后可有依据地收缩范围（这是本设计的取舍：宁可多列，不静默删除）。"
        )
    return diag


def explain_unimplemented(diag: Diagnosis, ruleset: RuleSet,
                          combine_enabled: bool = False,
                          search_scope: dict | None = None) -> list[str]:
    """说明哪些规则**没有**参与质量比对，避免用户误以为"没列出=不存在"。

    `combine_enabled` 决定口径文案：组合搜索开启时不能说"不搜索多事件叠加"，
    否则报告会自相矛盾（第 3b 节已经列了组合候选）。
    """
    no_mass = [r for r in ruleset.rules.values() if not (r.has_exact_mass or r.is_low_res_only)]
    msgs = [
        f"规则库共 {len(ruleset)} 条，其中 {len(ruleset.mass_matchable())} 条带可计算质量差参与本次比对；"
        f"{len(no_mass)} 条没有单一质量差（净式随结构/序列变化，或本身不是质量规则），**不能**按数值检索。"
    ]
    if no_mass:
        sample = sorted(no_mass, key=lambda r: r.rule_id)[:5]
        msgs.append(
            "例如：" + "、".join(f"{r.rule_id}({r.name_cn})" for r in sample)
            + "。这类规则需要按序列/条件查，见 `python -m dsh_core.cli --lookup`。"
        )
    if combine_enabled:
        msgs.append(
            "**搜索口径**：L1 比对「单条规则给出的单一质量差」（含其整数倍）；"
            "L2 在声明的事件预算内枚举**两个事件**的叠加（见第 3b 节与下面的范围统计）。"
            "3 个及以上事件、键编辑变换与中间体路径仍未搜索，"
            "因此「未列出」只在声明范围内成立，不能据此断言不存在。"
        )
    else:
        msgs.append(
            "**搜索口径**：本次只比对「单条规则给出的单一质量差」（含其整数倍），"
            "**未启用组合搜索**（`--no-combine`）。"
            "因此对 +31.99、+67.98 这类「可能由两个事件相加得到」的观测值，"
            "若库内没有对应的单条规则，本报告会显示无候选——这不等于该组合不存在。"
        )
    return msgs