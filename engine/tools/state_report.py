"""Human-facing Chinese report from the complete state-core result."""
from structure_state import encoded

def display_formula(formula):return ''.join(e+(str(n) if n!=1 else '') for e,n in sorted(formula.items()) if n)
def render_state_report(result):
    lines=['# 多肽副产物候选报告','',f"病例：{result['case'].get('case_id','未命名')}。质量吻合只提供候选，不能单独确证结构。",'',
           '## 输入与质量参照','',f"反应：{result['case']['reaction']['description']}",'']
    for r in result['reference_components']:
        lines.append(f"- 目标组分 {' + '.join(r['chain_ids'])}：{display_formula(r['formula'])}；元素质量 {r['element_mass']:.6f} Da；固定电荷 {r['intrinsic_charge']:+d}。")
    for r in result['starting_components']:
        lines.append(f"- 起始组分 {' + '.join(r['chain_ids'])}：{display_formula(r['formula'])}；元素质量 {r['element_mass']:.6f} Da。")
    for i,observation in enumerate(result['observations'],1):
        lines.extend(['',f'## 观测 {i} 与全部匹配候选','',f"原始观测：{encoded(observation['input'])}",'',
                      f"候选 {observation['candidate_count']} 个；条件冲突 {observation['conflict_count']} 个；离子假设 {observation['ion_hypothesis_count']} 个。",'',
                      '|候选|组分／分子式|理论观测值|误差 Da|条件|路径数|','|---|---|---:|---:|---|---:|'])
        if not observation['candidates']:lines.append('|—|搜索范围内暂无匹配|—|—|—|—|')
        for j,h in enumerate(observation['candidates'],1):
            condition='条件未知，保留' if h['condition']=='unknown' else '条件相容'
            lines.append(f"|{j}|{'+'.join(h['component']['chain_ids'])} / {display_formula(h['component']['formula'])}|{h['predicted']:.6f}|{h['error_da']:+.6f}|{condition}|{len(h['paths'])}|")
        for j,h in enumerate(observation['candidates'],1):
            lines.extend(['',f'### 候选 {j} 的结构与证据','',
                          f"质量口径：{h['mass_basis']}；离子假设：{h['ion_assumption']}；固定电荷：{h['component']['intrinsic_charge']:+d}。"])
            if h.get('propagated_error'):
                e=h['propagated_error'];lines.append(f"双峰差已传播误差：±{e['tolerance_da']:.6f} Da，采用 {e['error_model']} 模型。")
            for chain in h['component']['state']['chains']:
                p=chain['peptide'];lines.append(f"- 链 {chain['id']}：{'–'.join(p['sequence'])}；N端 {p['nterm']}；C端 {p['cterm']}；修饰 {encoded(p['modifications'])}。")
                annotations=[a for a in chain['annotations'] if a.get('kind')!='product_provenance']
                if annotations:lines.append('  连接／结构状态：'+encoded(annotations))
            if h['component']['state']['bonds']:lines.append('- 连接：'+encoded(h['component']['state']['bonds']))
            for k,path in enumerate(h['paths'][:3],1):
                if not path:lines.append(f'- 路径 {k}：目标／起始物本身，未假定附加反应。');continue
                lines.append(f"- 路径 {k}："+' → '.join(e['branch_id'] for e in path))
                for event in path:
                    lines.append(f"  - {event['rule_id']}，位点 {', '.join(event['sites'])}；{event['product']}；本步净式 {encoded(event['delta_formula'])}；来源 {encoded(event.get('source'))}；证据 {event['evidence_level']}。")
                    if event['missing_conditions']:lines.append('    待确认前提：'+'、'.join(event['missing_conditions'])+'。')
                    if event.get('elementary_exchange_steps'):lines.append('    桥网络仅列终点，具体交换步数未推定。')
            if len(h['paths'])>3:lines.append(f"其余 {len(h['paths'])-3} 条支持路径保留在完整JSON中。")
            lines.append('区分证据：结合合成步骤、保护历史、跨位点 MS/MS、同位素峰形、对照或标品核对；同式位点及连接不能只凭质量区分。')
            for note in h.get('source_notes',[]):
                lines.append(f"- {note['rule_id']}（{note['origin']}）的鉴别依据："+str(note['diagnostic_evidence']))
                lines.append('  同质量风险：'+str(note['competing_explanations']))
                lines.append('  工艺线索（应核对兼容性）：'+str(note['mitigation']))
        if observation['condition_conflicts']:
            lines.extend(['','### 已排除的质量匹配及理由',''])
            for h in observation['condition_conflicts']:
                reasons=sorted({x for p in h['paths'] for e in p for x in e['conflicting_conditions']})
                lines.append(f"- {display_formula(h['component']['formula'])}：明确否定的前提为 {'、'.join(reasons)}。")
        if observation.get('source_leads'):
            lines.extend(['','### 全书质量家族与附录线索','',
                          '以下是来源反查，不等于已执行的结构候选；参照、保护态或连接仍需核对。',
                          '|来源ID|质量家族／线索|净式／理论差|条件状态|来源|','|---|---|---|---|---|'])
            for lead in observation['source_leads']:
                name=str(lead['name']).replace('|','/').replace('\n',' ')
                mass=f"{lead['theoretical_delta']:+.6f}" if lead['theoretical_delta'] is not None else str(lead['nominal'])+'（原书整数线索）'
                condition='明确冲突：'+'、'.join(lead['conflicting_conditions']) if lead['context_status']=='context_conflict' else '来源条件支持，参照待核' if lead['context_status']=='context_supported' else '前提未知，参照待核'
                lines.append(f"|{lead['id']}|{name}|{mass}|{condition}|{encoded(lead['source']).replace('|','/')}|")
        if observation.get('dsh_source_review'):
            lines.extend(['','### DSH独立质量反查补充','',
                '以下仍是来源线索；使用同一规范组件和规则，条件以规范分支为准。不能把这些条目视为新执行结构。'])
            for lead in observation['dsh_source_review']:
                lines.append(f"- {lead['id']} {lead['name']}：{lead['theoretical_delta']:+.6f} Da，误差 {lead['error_da']:+.6f}；{lead['mass_grade']}；来源 {encoded(lead['source'])}。")
                lines.append('  区分依据：'+lead['diagnostic_evidence'])
    if result.get('fusion_audit'):
        a=result['fusion_audit'];lines.extend(['','## 融合与数据核对','',
            f"正式结构执行：{a['formal_executor']}；DSH反查 {a['dsh_rules']} 条历史来源；本体 {a['canonical_rules']} 条规则。",
            f"DSH独立产物分子式核算 {a['whole_product_formula_checks']} 项，通过：{a['whole_product_formula_checks_pass']}；数据差异：{len(a['dsh_data_problems'])}。"])
    lines.extend(['','## 搜索范围与能力边界',''])
    for s in result['search_scopes']:
        lines.append(f"- {s['origin']}：最多 {s['event_budget']} 次规则应用，生成 {s['nodes']} 个状态；耗时 {s['seconds']:.3f} 秒。")
        if s['truncated']:lines.append(f"  已截断：{s['truncation_reason']}；本次结果不能称完整。")
        if s['unsupported_transitions']:
            for gap in s['unsupported_transitions']:lines.append(f"  未执行 {gap['branch_id']}：{gap['reason']}。")
        if s['source_records_without_templates']:lines.append(f"  {len(s['source_records_without_templates'])} 条来源记录还没有结构执行模板，包括定性事实和仍需实现的模型；不以它们的质量线索冒充结构候选。")
    lines.extend(['',f"冻结知识包摘要：{result['library_hash']}。完整候选及全部路径保存在同次 JSON 结果中。",'',
                  '此报告采用反应家族与质量的实用粒度。未定义结构和复杂连接保留为来源线索；真实样品准确率未单独验收。'])
    return '\n'.join(lines)+'\n'
