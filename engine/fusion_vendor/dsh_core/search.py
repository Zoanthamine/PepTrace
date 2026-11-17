"""T06：**声明范围内完整**的组合搜索。

执行计划第 6 节要求：搜索阶段要先生成目标分子，再做单事件、序列片段与连接变化，
随后按当前分子状态组合；允许重复发生的不同位点要枚举；同一位点连续反应要有规则路径；
必须声明 event_budget、搜索层级、节点数、截断原因、未搜索范围；完整性 = 声明范围内完整。

本模块实现范围（与 search_scope.json 保持一致）：

- L0 目标分子：输入状态（序列 + 端基 + 修饰），质量由残基表独立计算。
- L1 单事件：规则库中所有**可数值检索**的规则（有精确净式/精确值，或只有整数 amu，
  或标注「及其倍数」），按适用性（残基前提）与条件四态判定后计算质量差。
- L2 两事件组合：任意两个 L1 事件叠加。组合的净式优先按元素相加，只有数值时才降级为
  数值相加（mass_only 语义）。同一位点冲突与重复上限按下面的保守策略判定。
- 未实现（明确写进 search_scope.json，不假装搜过）：3 个及以上事件、连接/环化/骨架重排的
  键编辑变换、多步中间体路径（中间体状态依赖的连续反应）、异构体/消旋的手性区分、
  片段缩合类的序列编辑。这些需要逐条的键编辑规则，规则库目前只有净式与文字描述。

设计取舍：
- 组合搜索**先不看观测值**地枚举可行事件，再统一比对质量；否则会漏掉"单事件都不匹配、
  两事件叠加才匹配"的情形（这正是组合搜索存在的理由）。
- 冲突判定是**保守**的：宁可多给一个并列候选（并标注冲突理由），也不静默删掉可能真实的
  组合。冲突只用于排除"同一端基被两次消耗"这类物理上不可能的情形。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

from .candidates import COND_LABEL, COND_RANK, _eval_conditions, _eval_residues
from .components import ComponentLibrary
from .formula import Formula, sum_formulas
from .matching import match_delta, ppm_error
from .rules import Rule, RuleSet
from .schema import CaseInput

ENGINE_VERSION = "dsh-t06-1.0.0"

# 搜索范围声明（同时写进 search_scope.json 与报告第 7 节）
UNSUPPORTED = [
    "3 个及以上事件的叠加（当前 event_budget=2）",
    "连接/环化/骨架重排的键编辑变换（规则库只给净式，未给可执行键编辑）",
    "多步中间体路径（中间体 → 后续反应的状态依赖链）",
    "异构体、消旋、二硫拓扑的手性/拓扑区分（净式相同处只作并列）",
    "多聚/二聚平衡关系（如 [2M−18] 写的是产物总质量），不当作单体质量差参与单事件或组合检索",
    "片段缩合、序列编辑（sequence_edit 恒为 null，规则库未给逐残基编辑指令）",
]


@dataclass(frozen=True)
class Event:
    """一个可参与组合的单事件（规则的一次应用）。"""

    rule_id: str
    name_cn: str
    unit_ratio: int
    site_key: str      # 'nterm' | 'cterm' | 'r:A' | 'unspecified'
    grade: str         # exact | nominal | exact_multiple | nominal_multiple
    delta_mass: float | None
    delta_formula: Formula | None


@dataclass
class Combination:
    """一条组合候选（1 个或 2 个事件）。"""

    events: tuple[Event, ...]
    predicted: float
    grade: str
    error_da: float
    error_ppm: float | None
    delta_formula: Formula | None
    condition_state: str
    condition_reason: str
    sites: list[dict] = field(default_factory=list)
    verdict: str = "retained"
    excluded_reason: str = ""
    singles_matched: tuple[bool, ...] = ()
    note: str = ""

    @property
    def key(self) -> str:
        return "+".join(f"{e.rule_id}x{e.unit_ratio}" if e.unit_ratio > 1 else e.rule_id
                        for e in self.events)

    def rank_key(self) -> tuple:
        grade = 0 if self.grade.startswith("exact") else 1
        err = abs(self.error_ppm) if self.error_ppm is not None else abs(self.error_da) * 1e3
        return (len(self.events), grade, -COND_RANK[self.condition_state], err)

    def as_row(self) -> dict:
        return {
            "combination_key": self.key,
            "events": [{"rule_id": e.rule_id, "name_cn": e.name_cn,
                        "unit_ratio": e.unit_ratio, "site": e.site_key} for e in self.events],
            "event_count": len(self.events),
            "predicted_delta_da": round(self.predicted, 6),
            "error_da": round(self.error_da, 6),
            "error_ppm": None if self.error_ppm is None else round(self.error_ppm, 2),
            "mass_grade": self.grade,
            "delta_formula": None if self.delta_formula is None else self.delta_formula.name(),
            "condition_state": self.condition_state,
            "condition_label": COND_LABEL[self.condition_state],
            "condition_reason": self.condition_reason,
            "sites": self.sites,
            "singles_matched": list(self.singles_matched),
            "verdict": self.verdict,
            "excluded_reason": self.excluded_reason,
            "note": self.note,
        }


@dataclass
class SingleHit:
    """一个与观测吻合的单事件（不进入 combinations，避免与 candidates 重复计数）。"""

    rule_id: str
    name_cn: str
    unit_ratio: int
    grade: str
    delta_mass: float
    error_da: float
    error_ppm: float | None
    condition_state: str
    condition_reason: str
    site_key: str

    @property
    def key(self) -> str:
        return f"{self.rule_id}x{self.unit_ratio}" if self.unit_ratio > 1 else self.rule_id

    def as_row(self) -> dict:
        return {
            "rule_id": self.rule_id,
            "name_cn": self.name_cn,
            "unit_ratio": self.unit_ratio,
            "predicted_delta_da": round(self.delta_mass, 6),
            "error_da": round(self.error_da, 6),
            "error_ppm": None if self.error_ppm is None else round(self.error_ppm, 2),
            "mass_grade": self.grade,
            "condition_state": self.condition_state,
            "condition_label": COND_LABEL[self.condition_state],
            "condition_reason": self.condition_reason,
            "site": self.site_key,
        }


@dataclass
class SearchResult:
    combinations: list[Combination] = field(default_factory=list)
    excluded: list[Combination] = field(default_factory=list)
    scope: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    # 单事件吻合的完整记录（不进 combinations，避免与 candidates 重复计数）。
    # 供独立穷举对照与审计使用：完整性比对必须是「单事件 ∪ 两事件」而不是只看组合。
    singles: list["SingleHit"] = field(default_factory=list)

    def retained_keys(self) -> list[str]:
        return [c.key for c in self.combinations]

    def all_hit_keys(self) -> list[str]:
        """单事件与两事件的全部吻合键（用于与独立穷举逐条对齐）。"""
        return [s.key for s in self.singles] + [c.key for c in self.combinations]


def _site_key(rule: Rule) -> str:
    sel = (rule.mass_status, rule.named_residues, rule.section)
    text = " ".join(str(x) for x in (rule.name_cn, rule.name_en, rule.conditions,
                                     rule.products, rule.susceptible))
    if rule.named_residues:
        return "r:" + ",".join(sorted(rule.named_residues))
    if "C 端" in text or "C端" in text:
        return "cterm"
    if "N 端" in text or "N端" in text:
        return "nterm"
    return "unspecified"


def _rule_text(rule: Rule) -> str:
    return " ".join(str(x) for x in (rule.name_cn, rule.name_en, rule.conditions,
                                     rule.products, rule.susceptible))


def _compare_state(a: str, b: str) -> tuple[str, str]:
    """两个事件的条件状态合成：取更差的一档（violated 最优先生效）。"""
    order = ["violated", "unknown", "weak", "satisfied"]
    worse = a if order.index(a) <= order.index(b) else b
    if a == b:
        return a, f"两个事件的条件判定均为 {COND_LABEL[a]}"
    return worse, f"两个事件条件判定分别为 {COND_LABEL[a]} / {COND_LABEL[b]}，按更保守的一档取 {COND_LABEL[worse]}"


def _conflict(rule_a: Rule, ea: int, rule_b: Rule, eb: int, case: CaseInput,
              lib: ComponentLibrary) -> str:
    """返回冲突理由；空串表示可以组合。保守策略：只排除物理上不可能的叠加。"""
    state = case.state
    ka, kb = _site_key(rule_a), _site_key(rule_b)

    # 同一端基被两条规则同时消耗 → 不可能
    for key in ("nterm", "cterm"):
        if ka == key and kb == key:
            return f"两条规则都要消耗同一个 {key} 位点，同一位点连续反应需要规则路径，当前引擎不建模"

    # 同一条规则重复应用：只有明确允许重复、且肽里有多个可用位点时才允许
    if rule_a.rule_id == rule_b.rule_id:
        if (rule_a.nominal_multiples or 1) <= 1:
            return "同一规则的重复应用上限为 1（书本条标的是单价反应），当前不支持同一位点连续反应"
        residues = [x for x in rule_a.named_residues if x != "pGlu"]
        if not residues:
            return ""
        if all(len(state.positions_of(x)) >= 2 for x in residues):
            return ""      # 多个不同位点，允许各一次
        return ("该规则涉及残基在序列中只有一处，两次应用会落在同一位点；"
                "同一位点连续反应需要规则路径，当前不支持简单相加")

    # 同一残基类型、肽里只有一处 → 两条规则会争同一个位点
    if ka.startswith("r:") and kb.startswith("r:"):
        shared = set(ka[2:].split(",")) & set(kb[2:].split(","))
        if shared and ea == eb == 1:
            single = [x for x in shared if len(state.positions_of(x)) == 1]
            if single:
                return (f"两条规则都要作用于同一处 {'/'.join(sorted(single))}（序列中只有一处），"
                        "同一位点连续反应需要规则路径，简单相加会给出不存在的组合")
    return ""


def _event_window(ev: "Event", tol_da: float) -> float:
    """单个事件的匹配窗口：tol·k，整数口径再叠加 0.5·k 的舍入箱。"""
    k = ev.unit_ratio
    return tol_da * k + (0.5 * k if ev.grade.startswith("nominal") else 0.0)


def _fits(predicted: float, observed: float, tol_da: float,
          grade: str, unit_ratio: int = 1) -> tuple[bool, float, float | None]:
    """按档位判定是否落在容差内，窗口规则与 matching.match_delta 保持一致。

    exact  : |实测−理论| ≤ tol
    nominal: |实测−理论| ≤ 0.5 + tol（整数 amu 口径与精确口径绝不混同）
    multiple: 窗口按倍数放缩（k·tol 或 0.5k + k·tol），与 match_delta 同规则。
    """
    err = observed - predicted
    if grade.startswith("exact"):
        window = tol_da * unit_ratio
    else:
        window = 0.5 * unit_ratio + tol_da * unit_ratio
    if abs(err) <= window:
        return True, err, ppm_error(observed, predicted)
    return False, err, None


def _event_delta(rule: Rule, unit_ratio: int) -> tuple[Formula | None, float | None, str]:
    """返回 (净式, 质量, 档位)。净式优先，只有整数 amu 时降级 nominal。"""
    if rule.delta is not None:
        f = rule.delta * unit_ratio
        grade = "exact_multiple" if unit_ratio > 1 else "exact"
        return f, f.mass(), grade
    if rule.nominal is not None:
        return None, rule.nominal * unit_ratio, ("nominal_multiple" if unit_ratio > 1 else "nominal")
    return None, None, "none"


def search_combinations(case: CaseInput, ruleset: RuleSet, lib: ComponentLibrary,
                        max_events: int = 2, max_nodes: int = 200_000,
                        max_results: int = 5_000) -> SearchResult:
    """在声明范围内穷举 1–2 个事件的组合，返回质量吻合的组合与范围统计。"""
    res = SearchResult()
    obs = case.observed
    target_formula = case.state.formula(lib)
    target_mass = target_formula.mass()

    if obs.kind == "delta":
        observed = float(obs.value)
        delta_source = "输入即中性质量差"
    else:
        # 绝对质量（neutral_mass / mz / mz_pair 归一后）先反推成"产物−起始物"
        neutral = obs.neutral if obs.neutral is not None else float(obs.value)
        observed = float(neutral) - target_mass
        delta_source = (f"由绝对质量 {neutral:.6f} Da 减去起始物质量 {target_mass:.6f} Da 反推"
                        f"（等价于 ΔM={observed:+.6f} Da；绝对质量路径本身由匹配层按产物质量计算）")
    tol = obs.effective_tol_da

    # ---- 1) 生成 L1 候选事件（不看观测值，先按适用性与条件筛）----
    events: list[Event] = []
    applicability: dict[str, str] = {}
    condition_of: dict[str, tuple[str, str]] = {}
    skipped_no_mass = 0
    skipped_residue = 0
    skipped_polymer: list[str] = []
    skipped_zero_mass: list[str] = []
    rule_list = sorted(ruleset.rules.values(), key=lambda r: r.rule_id)
    for rule in rule_list:
        if rule.polymer_relation:
            # [2M−18] 这类是「产物总质量」的多聚关系，不是单体质量差，不能当单事件叠加
            skipped_polymer.append(f"{rule.rule_id}（{rule.polymer_relation}）")
            continue
        if not (rule.has_exact_mass or rule.is_low_res_only):
            skipped_no_mass += 1
            continue
        app, reason, sites = _eval_residues(rule, case.state)
        if app != "applicable":
            skipped_residue += 1
            continue
        cond_state, cond_reason = _eval_conditions(rule, case)
        if cond_state == "violated":
            res.excluded.append(Combination(
                events=(Event(rule.rule_id, rule.name_cn, 1, _site_key(rule), "excluded",
                              rule.delta_mass, rule.delta),),
                predicted=rule.delta_mass or 0.0, grade="excluded", error_da=0.0,
                error_ppm=None, delta_formula=rule.delta, condition_state=cond_state,
                condition_reason=cond_reason, sites=sites, verdict="excluded",
                excluded_reason=cond_reason))
            continue
        applicability[rule.rule_id] = app
        condition_of[rule.rule_id] = (cond_state, cond_reason)
        ratios = [1]
        if rule.nominal_multiples:
            ratios = list(range(1, min(int(rule.nominal_multiples), 8) + 1))
        for ratio in ratios:
            formula, mass, grade = _event_delta(rule, ratio)
            if mass is None:
                # 净质量差 0 的异构类规则（delta_formula 与 nominal 都为空、delta_mass=0）：
                # 它们不改变质量，进组合层只会把每个质量解膨胀数十倍且无法用质量区分，
                # 故不参与质量组合；单事件层仍以 ΔM=0 并列呈现（isomeric_no_mass_change）。
                if rule.delta_mass == 0.0:
                    skipped_zero_mass.append(rule.rule_id)
                continue
            events.append(Event(rule.rule_id, rule.name_cn, ratio, _site_key(rule),
                                grade, mass, formula))

    # ---- 2) 单事件 ----
    nodes = 0
    truncated = False
    truncation_reasons: list[str] = []
    single_matched: dict[tuple[str, int], bool] = {}
    single_hits: list[tuple[Event, float, float]] = []
    for ev in events:
        nodes += 1
        if nodes > max_nodes:
            truncated = True
            truncation_reasons.append(f"节点预算 {max_nodes} 用尽：未搜索全部单事件与组合")
            break
        hit, err_da, err_ppm = _fits(ev.delta_mass, observed, tol, ev.grade, ev.unit_ratio)
        single_matched[(ev.rule_id, ev.unit_ratio)] = hit
        if not hit:
            continue
        # 单事件吻合已由 diagnose() 的 candidates 完整呈现；combinations 只放 ≥2 事件，
        # 避免同一结果在两处重复计数（曾因此让 §3b 表里出现 event_count=1 的行）。
        single_hits.append((ev, err_da, err_ppm))

    # ---- 3) 两事件组合 ----
    pair_examined = 0
    pair_conflict = 0
    conflict_keys: list[str] = []
    if max_events >= 2 and not truncated:
        rule_by_id = ruleset.rules
        for ev_a, ev_b in combinations(events, 2):
            if ev_a.rule_id == ev_b.rule_id and ev_a.unit_ratio == ev_b.unit_ratio:
                continue
            # 规则内顺序：同一条规则的不同倍数不当作两个事件
            if ev_a.rule_id == ev_b.rule_id:
                continue
            pair_examined += 1
            nodes += 1
            if nodes > max_nodes:
                truncated = True
                truncation_reasons.append(f"节点预算 {max_nodes} 用尽：两事件组合未穷举完")
                break
            conflict = _conflict(rule_by_id[ev_a.rule_id], ev_a.unit_ratio,
                                 rule_by_id[ev_b.rule_id], ev_b.unit_ratio, case, lib)
            if conflict:
                pair_conflict += 1
                if len(conflict_keys) < 20_000:
                    conflict_keys.append(f"{ev_a.rule_id}x{ev_a.unit_ratio}+{ev_b.rule_id}x{ev_b.unit_ratio}"
                                         if (ev_a.unit_ratio > 1 or ev_b.unit_ratio > 1)
                                         else f"{ev_a.rule_id}+{ev_b.rule_id}")
                continue
            if ev_a.delta_formula is not None and ev_b.delta_formula is not None:
                total_formula = sum_formulas([ev_a.delta_formula, ev_b.delta_formula])
                total_mass = total_formula.mass()
                grade = ("exact" if ev_a.grade.startswith("exact") and ev_b.grade.startswith("exact")
                         else "nominal")
            else:
                total_formula = None
                total_mass = ev_a.delta_mass + ev_b.delta_mass
                grade = "nominal"
            # 组合窗口 = 两个单事件窗口之和（误差合成，不取平均）
            window = _event_window(ev_a, tol) + _event_window(ev_b, tol)
            err = observed - total_mass
            if abs(err) > window:
                continue
            err_ppm = ppm_error(observed, total_mass)
            st_a, rs_a = condition_of[ev_a.rule_id]
            st_b, rs_b = condition_of[ev_b.rule_id]
            cstate, creason = _compare_state(st_a, st_b)
            singles = (single_matched.get((ev_a.rule_id, ev_a.unit_ratio), False),
                       single_matched.get((ev_b.rule_id, ev_b.unit_ratio), False))
            note = ("组合事件：单事件都不匹配，只有叠加后才吻合" if not any(singles)
                    else "组合事件（其中至少一个单事件也吻合，组合是另一种解释）")
            if grade == "nominal" and abs(err) > 0.05:
                # 整数口径规则（书里只给 +N amu）自带 ±0.5 Da 舍入箱；当观测是高分辨时，
                # 这类命中的误差可能远大于高分辨精度，只能当粗筛线索，不能当确认证据。
                note += "；整数口径放宽命中，误差大于高分辨精度，仅作粗筛线索，需 MS/MS 或标品确认"
            rule_a = rule_by_id[ev_a.rule_id]
            _, _, sites_a = _eval_residues(rule_a, case.state)
            rule_b = rule_by_id[ev_b.rule_id]
            _, _, sites_b = _eval_residues(rule_b, case.state)
            res.combinations.append(Combination(
                events=(ev_a, ev_b), predicted=total_mass, grade=grade,
                error_da=err, error_ppm=err_ppm, delta_formula=total_formula,
                condition_state=cstate, condition_reason=creason,
                sites=sites_a + sites_b, singles_matched=singles, note=note))

    # ---- 4) 结果排序与上限 ----
    res.combinations.sort(key=lambda c: (c.rank_key(), c.key))
    truncated_by_results = False
    if len(res.combinations) > max_results:
        truncated_by_results = True
        truncation_reasons.append(f"结果数上限 {max_results}：超出部分未列出（已排序，截断保留最优）")
        res.combinations = res.combinations[:max_results]

    product_mass = target_mass + observed
    res.scope = {
        "engine": ENGINE_VERSION,
        "enabled": True,
        "max_events": max_events,
        "event_budget": {
            "max_events": max_events,
            "max_nodes": max_nodes,
            "max_results": max_results,
        },
        "levels": ["L0 目标分子", "L1 单事件", "L2 两事件组合"],
        "counts": {
            "rules_total": len(rule_list),
            "rules_mass_matchable": len(rule_list) - skipped_no_mass,
            "rules_skipped_no_mass": skipped_no_mass,
            "rules_skipped_residue": skipped_residue,
            "events_generated": len(events),
            "nodes_expanded": nodes,
            "singles_matched": sum(1 for v in single_matched.values() if v),
            "pairs_examined": pair_examined,
            "rules_skipped_polymer_relation": len(skipped_polymer),
            "rules_skipped_zero_mass_isomeric": len(skipped_zero_mass),
            "zero_mass_isomeric_rules": skipped_zero_mass[:40],
            "zero_mass_isomeric_rules_truncated": len(skipped_zero_mass) > 40,
            "polymer_relation_rules": skipped_polymer,
            "pairs_rejected_by_site_conflict": pair_conflict,
            "conflict_keys_sample": conflict_keys,
            "conflict_keys_truncated": pair_conflict > len(conflict_keys),
            "combinations_retained": len(res.combinations),
            "combinations_excluded": len(res.excluded),
        },
        "observed_delta_da": round(observed, 6),
        "delta_source": delta_source,
        "tolerance_da": round(tol, 6),
        "target_mass": round(target_mass, 6),
        "implied_product_mass": round(product_mass, 6),
        "truncated": truncated or truncated_by_results,
        "truncation_reasons": truncation_reasons,
        "unsupported": UNSUPPORTED,
        "completeness_note": (
            "本结果的完整性只覆盖 L1 单事件与 L2 两事件组合；触顶时 truncated=true 且列出未搜索范围。"
            "净质量差为 0 的异构类规则不参与质量组合（不改变质量，只会让每个质量解膨胀；"
            "它们在单事件层以 ΔM=0 并列显示）；3 个以上事件、键编辑与中间体路径未搜索（见 unsupported），"
            "不能据此断言“不存在其他副反应”。"
        ),
    }
    if truncated:
        res.notes.append("搜索触顶被截断：候选集合不完整，报告须显式提示。")
    # 单事件吻合记录（供审计/穷举对照；不进 §3b 组合表，避免与 candidates 重复计数）
    res.singles = [
        SingleHit(
            rule_id=ev.rule_id, name_cn=ev.name_cn, unit_ratio=ev.unit_ratio,
            grade=ev.grade, delta_mass=float(ev.delta_mass),
            error_da=err_da, error_ppm=err_ppm,
            condition_state=condition_of.get(ev.rule_id, ("unknown", ""))[0],
            condition_reason=condition_of.get(ev.rule_id, ("unknown", ""))[1],
            site_key=ev.site_key,
        )
        for ev, err_da, err_ppm in single_hits
    ]
    return res