"""中文报告渲染。

报告顺序刻意按"诊断者的思考顺序"排：
先确认我理解了什么（序列/端基/质量）→ 观测口径 → 候选（含位点与条件）→
同式并列 → 反证 → 我**没**能比对的部分 → 下一步验证。

诚实性约束
- 书中观察（origin=book）与编辑补充（origin=supplement）必须分开标注，补充条目不得被当成原书数据引用。
- 未知就写"未知"，不编条件、不编位点。
"""
from __future__ import annotations

from .candidates import COND_LABEL, Diagnosis, explain_unimplemented
from .rules import RuleSet

ORIGIN_LABEL = {"book": "书中", "supplement": "补充"}


def _fmt_pps(pps: list[int]) -> str:
    if not pps:
        return "—"
    if len(pps) == 1:
        return f"p.{pps[0]}"
    return f"p.{pps[0]}–{pps[-1]}" if pps == list(range(pps[0], pps[-1] + 1)) else "p." + ",".join(str(p) for p in pps)


def _sites_text(sites: list[dict]) -> str:
    if not sites:
        return "位点未定（规则未限定残基，或需按序列人工判断）"
    bits = []
    for s in sites:
        pos = "/".join(str(p) for p in s["positions"]) if s["positions"] else "?"
        bits.append(f"{s['residue']}@{pos}")
    return "、".join(bits)


def render_markdown(diag: Diagnosis, ruleset: RuleSet, top_k: int | None = None,
                    search=None) -> str:
    case = diag.case
    state = case.state
    lines: list[str] = []
    add = lines.append

    add(f"# 多肽副产物诊断报告：{case.case_id}")
    add("")
    add("> 本报告由 DSH 版诊断引擎（`vnext_dsh`）生成。候选来自 Yang 2016《Side Reactions in Peptide Synthesis》"
        "全书 382 条规则 + 63 条编辑补充；每条都带规则号与印刷页码，可回查。")
    add("")

    # ---------------- 1 输入理解 ----------------
    add("## 1. 输入理解（先确认我算的是不是你要的东西）")
    add("")
    add(f"- 规范态：`{state.describe()}`")
    add(f"- 序列长度：{state.length}；残基组成：{''.join(sorted(set(state.sequence)))}")
    add(f"- 起始物中性单同位素质量：**{diag.target_mass:.6f} Da**（净式 {diag.target_formula_name}）")
    if state.folded:
        add(f"- 已折叠的等价写法：{'；'.join(state.folded)}")
    for w in case.warnings:
        add(f"- ⚠️ {w['message']}")
    add("")

    # ---------------- 2 观测 ----------------
    add("## 2. 观测口径")
    add("")
    add(f"- 输入口径：{diag.observed_desc}")
    if diag.observed_delta is not None:
        add(f"- 归一到中性质量差：**ΔM = {diag.observed_delta:+.6f} Da**（= 产物 − 起始物，中性）")
        add(f"- 检索容差：±{diag.tol_da:.4f} Da")
    obs = case.observed
    if obs.pair_derivation:
        add(f"- {obs.pair_derivation}")
        add("- 两峰求差已核对电荷与加合：电荷不同按各自 z 换算中性质量后再相减，误差按方和根传播。")
    if case.fixed_charge:
        add(f"- 固定电荷：{case.fixed_charge:+d}（永久带电基团）；m/z ↔ 中性质量换算**不加质子**，"
            f"与 [M+{abs(case.fixed_charge)}H] 口径相差约 {1.007276466621 * abs(case.fixed_charge) / max(1, abs(obs.z or 1)):.4f} Da/z。")
    add(f"- 质量差符号约定：**ΔM = M(产物) − M(起始物)**；m/z 差必须乘电荷后再比较。")
    add("")

    # ---------------- 3 候选 ----------------
    shown = diag.candidates if top_k is None else diag.candidates[:top_k]
    add(f"## 3. 质量吻合的候选（共 {len(diag.candidates)} 条，此处列 {len(shown)} 条）")
    add("")
    if not shown:
        add("_本次没有保留任何质量吻合的候选。_" if not diag.excluded
            else "_质量命中的候选全部被排除，理由见第 5 节。_")
        add("")
    else:
        add("| # | 规则 | 结构族/反应 | 预测 ΔM | 误差 | 质量档 | 条件判定 | 候选位点 | 出处 |")
        add("|---|---|---|---|---|---|---|---|---|")
        for i, c in enumerate(shown, 1):
            err = f"{c.match.error_da:+.4f} Da"
            if c.match.error_ppm is not None:
                err += f" ({c.match.error_ppm:+.0f} ppm)"
            add(
                f"| {i} | {c.rule_id} | {c.rule.name_cn} | {c.match.predicted:+.6f} | {err} | "
                f"{c.match.grade} | {COND_LABEL[c.condition_state]} | {_sites_text(c.sites)} | "
                f"{ORIGIN_LABEL.get(c.rule.origin, c.rule.origin)}·{_fmt_pps(list(c.rule.print_pages))} |"
            )
        add("")
        add("### 候选明细（含判定依据与冲突）")
        add("")
        for c in shown:
            r = c.rule
            add(f"#### {c.rule_id} {r.name_cn}")
            add("")
            add(f"- 预测：ΔM = {c.match.predicted:+.6f} Da"
                + (f"（净式 {r.delta.name()}）" if r.delta else "")
                + (f"；预计产物质量 {c.expected_product_mass:.6f} Da" if c.expected_product_mass else ""))
            add(f"- 误差：实测 − 理论 = {c.match.error_da:+.6f} Da"
                + (f" = {c.match.error_ppm:+.1f} ppm" if c.match.error_ppm is not None else ""))
            if c.product_formula is not None:
                add(f"- 预计产物净式：{c.product_formula.name()}（质量 {c.product_formula.mass():.6f} Da）")
            if c.expected_mz is not None:
                add(f"- 预计 m/z：**{c.expected_mz:.6f}**（z={obs.z}，加合 {obs.adduct}"
                    + (f"，固定电荷 {case.fixed_charge:+d}" if case.fixed_charge else "") + "）")
            if c.isotope_lines:
                iso = "；".join(
                    f"M+{row['offset']} {row['mz']:.4f}（相对强度 {row['rel_intensity']:.2f}）"
                    for row in c.isotope_lines)
                add(f"- 同位素包络（用于核对实测同位素间距）：{iso}")
            add(f"- 质量档：{c._match_note()}；库内质量状态 `{r.mass_status}`")
            add(f"- 条件：{c.condition_reason}")
            add(f"- 序列适用性：{c.applicability_reason}")
            add(f"- 候选位点：{_sites_text(c.sites)}")
            if r.mass_expression:
                add(f"- 书中质量表达：{r.mass_expression}")
            if r.diagnostic_evidence:
                add(f"- 诊断证据：{r.diagnostic_evidence}")
            if r.competing:
                add(f"- 可能的竞争解释/同式风险：{r.competing}")
            if r.mitigation:
                add(f"- 缓解/验证：{r.mitigation}")
            if r.section:
                add(f"- 出处：第 {r.chapter} 章 {r.section}，{_fmt_pps(list(r.print_pages))}"
                    + (f"，图 {', '.join(r.figures)}" if r.figures else ""))
            add("")

    # ---------------- 3b 组合事件（T06） ----------------
    if search is not None:
        add("### 3b. 组合事件（两事件叠加；单事件口径之外的解释）")
        add("")
        if not search.combinations and not search.excluded:
            add("_声明范围内没有两事件组合落在容差内。_")
            add("")
        else:
            combos = [c for c in search.combinations if len(c.events) > 1]
            if not combos:
                add("_两事件组合在声明范围内没有质量吻合项；下面只有单事件项被组合层复述，已并入第 3 节。_")
            else:
                sample = combos if top_k is None else combos[:top_k]
                add(f"两事件组合质量吻合 {len(combos)} 条，此处列 {len(sample)} 条。"
                    "「单事件吻合」列为 False 表示两个单独事件都不在容差内、只有叠加后才吻合——"
                    "这类解释在单事件搜索引擎里是查不到的。")
                add("")
                add("| # | 组合 | 预测 ΔM | 误差 | 质量档 | 条件 | 单事件吻合 | 说明 |")
                add("|---|---|---|---|---|---|---|---|")
                for i, c in enumerate(sample, 1):
                    err = f"{c.error_da:+.4f} Da"
                    if c.error_ppm is not None:
                        err += f" ({c.error_ppm:+.0f} ppm)"
                    flags = "/".join("是" if x else "否" for x in c.singles_matched)
                    add(f"| {i} | {c.key} | {c.predicted:+.6f} | {err} | {c.grade} | "
                        f"{COND_LABEL[c.condition_state]} | {flags} | {c.note} |")
                add("")
            if search.excluded:
                add(f"另有 {len(search.excluded)} 条单事件因**用户明确否定的条件**被组合层排除（并入第 5 节口径）。")
                add("")

    # ---------------- 4 同式并列 ----------------
    add("## 4. 同净式/同质量候选（必须并列，不能只留一个）")
    add("")
    if not diag.groups:
        add("_本次保留的候选之间没有相同净式/相同质量的情况。_")
        add("")
    else:
        add("| 净式 | ΔM | 相同质量的规则 |")
        add("|---|---|---|")
        for g in diag.groups:
            add(f"| {g['delta_formula']} | {g['delta_mass']:+.6f} | "
                + "；".join(f"{rid} {nm}" for rid, nm in zip(g["rule_ids"], g["names"])) + " |")
        add("")
        add("> 质量相同≠结构相同。报告必须并列，再用保留时间、二级谱、加合模式、同位素峰等区分。")
        add("")

    # ---------------- 5 已排除 ----------------
    if not diag.excluded:
        add("## 5. 已排除的候选")
        add("")
        add("_本次没有被排除的质量吻合候选。_")
        add("")
    else:
        add(f"## 5. 已排除的候选（共 {len(diag.excluded)} 条，附理由）")
        add("")
        add("| 规则 | 反应 | 预测 ΔM | 误差 | 排除依据 |")
        add("|---|---|---|---|---|")
        for c in diag.excluded[: top_k or 30]:
            add(f"| {c.rule_id} | {c.rule.name_cn} | {c.match.predicted:+.6f} | {c.match.error_da:+.4f} Da | {c.excluded_reason} |")
        add("")

    # ---------------- 6 反证 ----------------
    if not diag.counter_evidence:
        add("## 6. 条件吻合但质量不符（反证）")
        add("")
        add("_本次没有需要单列的反证候选。_")
        add("")
    else:
        add("## 6. 条件吻合但质量不符（反证，用于否定「看起来像书里某条」的猜测）")
        add("")
        add("| 规则 | 反应 | 预测 ΔM | 与观测差距 | 说明 |")
        add("|---|---|---|---|---|")
        for x in diag.counter_evidence:
            add(f"| {x['rule_id']} | {x['name_cn']} | {x['predicted_delta_da']:+.6f} | {x['gap_da']:+.4f} Da | {x['why']} |")
        add("")

    # ---------------- 7 覆盖边界 ----------------
    add("## 7. 本次比对的边界（没列出的不等于不存在）")
    add("")
    for msg in explain_unimplemented(diag, ruleset,
                                     combine_enabled=bool(search is not None and search.scope.get("enabled")),
                                     search_scope=None if search is None else search.scope):
        add(f"- {msg}")
    if diag.data_problems:
        add(f"- ⚠️ 知识库自检发现 {len(diag.data_problems)} 处数据不一致（净式重算与记录值不符），已单列，不参与排序：")
        for p in diag.data_problems[:8]:
            add(f"  - `{p['rule_id']}`：{p['detail']}")
        if len(diag.data_problems) > 8:
            add(f"  - 其余 {len(diag.data_problems) - 8} 处见 JSON 输出的 `data_problems`。")
    add(f"- 统计：规则 {diag.stats['rule_total']} 条；可数值检索 {diag.stats['mass_matchable']} 条；"
        f"本次质量命中 {diag.stats['mass_matched']} 条 → 保留 {diag.stats['retained']}，排除 {diag.stats['excluded']}"
        f"（其中条件冲突 {diag.stats['excluded_condition_conflict']}，残基不符 {diag.stats['excluded_residue_missing']}）。")
    if search is not None:
        sc = search.scope
        if sc.get("enabled"):
            cnt = sc.get("counts", {})
            add(f"- 组合搜索范围（T06，引擎 `{sc.get('engine')}`）：搜索层级 {' → '.join(sc.get('levels', []))}；"
                f"生成事件 {cnt.get('events_generated')} 个，展开节点 {cnt.get('nodes_expanded')} 个"
                f"（上限 {sc.get('event_budget', {}).get('max_nodes')}），"
                f"两事件组合检查 {cnt.get('pairs_examined')} 对，位点冲突剔除 {cnt.get('pairs_rejected_by_site_conflict')} 对。")
            if sc.get("truncated"):
                add(f"- ⚠️ **组合搜索被截断**：{'；'.join(sc.get('truncation_reasons') or ['原因未记录'])}。"
                    "截断时候选集合不完整，不能据此断言「不存在其他副反应」。")
            else:
                add("- 组合搜索**未触顶**（在声明的 event_budget 内穷举完成）。")
            add(f"- 未实现/未搜索范围：{'；'.join(sc.get('unsupported', []))}")
            add(f"- 完整性口径：{sc.get('completeness_note')}")
        else:
            add(f"- 组合搜索未启动：{'; '.join(search.notes) or sc.get('reason')}")
    add("")

    # ---------------- 8 结论 ----------------
    add("## 8. 结论与下一步")
    add("")
    if diag.candidates:
        top = diag.candidates[0]
        add(f"- 最优先核对：**{top.rule_id} {top.rule.name_cn}**"
            f"（{COND_LABEL[top.condition_state]}，误差 {top.match.error_da:+.4f} Da）。")
        sat = [c for c in diag.candidates if c.condition_state == "satisfied"]
        unk = [c for c in diag.candidates if c.condition_state == "unknown"]
        if sat:
            add(f"- 有 {len(sat)} 条候选的条件已被你报告的信息**支持**，优先看这些。")
        if unk:
            add(f"- 有 {len(unk)} 条只是「质量吻合 + 条件未知」。要收缩范围，请补：所用试剂（如 TFA/EDT/TIS/哌啶）、"
                "脱保护/裂解步骤、温度与时长、是否含清除剂；以及是否**明确没有**用某个试剂。")
    elif diag.excluded:
        add(f"- 质量吻合的 {len(diag.excluded)} 条候选**全部被排除**（理由见第 5 节）。"
            "若你确信其中一个不该被排除，说明输入里的条件或序列信息需要更正——"
            "这正是本引擎要求「阴性条件必须显式报出」的原因。")
    else:
        add("- 没有质量吻合的候选，先核对输入（端基、修饰、电荷）与观测峰归属。")
    add("")
    add("- 验证手段（按可操作性排序）：同批次做**加标**实验、换电荷态看质量是否一致、"
        "看同位素峰形（是否引入 S/Br/Cl）、二级谱找特征 b/y 离子、必要时合成标准品共流出。")
    add("")
    for n in diag.notes:
        add(f"> {n}")
        add("")
    add("---")
    add("")
    add("说明：本报告是**可追溯的诊断辅助**，不是结构鉴定结论。书中数值、编辑补充数值、"
        "以及本引擎按净式重算的数值三者分开标注；结构图未经视觉核验的部分已在规则库中标为 `structure_required`。")
    return "\n".join(lines)


def render_json(diag: Diagnosis, ruleset: RuleSet, top_k: int | None = None,
                search=None) -> dict:
    shown = diag.candidates if top_k is None else diag.candidates[:top_k]
    case = diag.case
    return {
        "engine": "vnext_dsh",
        "case_id": case.case_id,
        "input": {
            "normalized_state": case.state.describe(),
            "sequence": case.state.sequence,
            "nterm_group": case.state.nterm_group,
            "cterm_amide": case.state.cterm_amide,
            "modifications": [
                {"site": m.site, "group": m.group} for m in case.state.modifications
            ],
            "folded_equivalences": list(case.state.folded),
            "conditions": {
                "reported": sorted(case.conditions.reported),
                "denied": sorted(case.conditions.denied),
                "tags": sorted(case.conditions.tags),
                "stage": case.conditions.stage,
            },
            "observed": {
                "kind": case.observed.kind,
                "value": case.observed.value,
                "z": case.observed.z,
                "adduct": case.observed.adduct,
                "tol_da": case.observed.tol_da,
                "tol_ppm": case.observed.tol_ppm,
                "normalized_neutral_mass": case.observed.neutral,
            },
            "warnings": case.warnings,
        },
        "target": {
            "monoisotopic_neutral_mass_da": round(diag.target_mass, 6),
            "formula": diag.target_formula_name,
        },
        "observed_delta_da": None if diag.observed_delta is None else round(diag.observed_delta, 6),
        "tolerance_da": round(diag.tol_da, 6),
        "candidates": [c.as_row() for c in shown],
        "excluded": [c.as_row() for c in diag.excluded],
        "same_delta_groups": diag.groups,
        "counter_evidence": diag.counter_evidence,
        "coverage_note": explain_unimplemented(
            diag, ruleset,
            combine_enabled=bool(search is not None and search.scope.get("enabled")),
            search_scope=None if search is None else search.scope),
        "combinations": [] if search is None else [c.as_row() for c in search.combinations],
        "combinations_excluded": [] if search is None else [c.as_row() for c in search.excluded],
        # 单事件吻合记录：与 candidates 同源，单独列出便于审计「单事件 ∪ 两事件」的完整性
        "search_singles": [] if search is None else [s.as_row() for s in search.singles],
        "search_scope": None if search is None else search.scope,
        "search_notes": [] if search is None else list(search.notes),
        "data_problems": diag.data_problems,
        "stats": diag.stats,
        "notes": diag.notes,
        "kb": {
            "rule_total": len(ruleset),
            "source": "Yang 2016, Side Reactions in Peptide Synthesis (全书 14 章) + 编辑补充条目",
            "generated_by": "vnext_dsh",
        },
    }