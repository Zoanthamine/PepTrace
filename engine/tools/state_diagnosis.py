"""Codex state-search CLI: neutral/mz/delta/paired observations with complete candidate files."""
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from case_input import check_keys, string_set, add_starting_group_facts
from carbonyl_components import normalize_carbonyl_donors
from knowledge_source import load_runtime
from observation_model import finite,formula_mass,hypotheses,peak_pair
from state_search import enumerate_states
from structure_state import normalize_molecule,describe,encoded,physical_fingerprint,fingerprint
from source_leads import mass_source_leads

def prepare_context(case,target,starting):
    reaction=case['reaction'];context=case.get('context',{})
    check_keys(reaction,{'description','stage_tags','stage_tags_absent','acyl_donors','carbonyl_donors'},'reaction')
    check_keys(context,{'features_present','features_absent','stage_tags_absent'},'context')
    if not isinstance(reaction.get('description'),str) or not reaction['description'].strip():raise ValueError('Provide the actual reaction description')
    present={'stage_tags':string_set(reaction.get('stage_tags',[]),'stage_tags'),
             'features':string_set(context.get('features_present',[]),'features'),'residues':set()}
    absent={'stage_tags':string_set(reaction.get('stage_tags_absent',[]),'stage_tags_absent')|string_set(context.get('stage_tags_absent',[]),'stage_tags_absent'),
            'features':string_set(context.get('features_absent',[]),'features_absent'),'residues':set()}
    if present['stage_tags']&absent['stage_tags'] or present['features']&absent['features']:raise ValueError('Contradictory context facts')
    from case_input import normalize_peptide
    present['acyl_donors']=[normalize_peptide(x) for x in reaction.get('acyl_donors',[])]
    if any(x['cterm']!='acid' for x in present['acyl_donors']):raise ValueError('Acyl donor requires an explicit free carboxylic acid')
    present['carbonyl_donors']=normalize_carbonyl_donors(reaction.get('carbonyl_donors',[]))
    if starting and len(starting['chains'])==len(target['chains'])==1:
        add_starting_group_facts(starting['chains'][0]['peptide'],target['chains'][0]['peptide'],(present,absent))
    return present,absent

def reference_components(case,observation,target,starting):
    ref=observation.get('reference','intended_product')
    if ref not in {'intended_product','starting_material'}:raise ValueError('Unknown observation reference')
    selected=target if ref=='intended_product' else starting
    if selected is None:raise ValueError('Starting-material comparison requires the starting structure')
    result=describe(selected)
    ids=observation.get('reference_chain_ids')
    if ids is not None:
        if not isinstance(ids,list) or not ids:raise ValueError('reference_chain_ids requires a nonempty chain ID list')
        result=[x for x in result if set(x['chain_ids'])==set(ids)]
        if not result:raise ValueError('Declared reference chains do not form one connected component')
    return result

def observation_options(observation):
    kind=observation.get('kind')
    if kind in {'delta_mono','delta_approx'}:
        check_keys(observation,{'kind','value','tolerance','tolerance_unit','reference','reference_chain_ids','mass_basis'},'relative observation')
        value=finite(observation.get('value'),'delta value');tolerance=finite(observation.get('tolerance'),'delta tolerance',True)
        if observation.get('tolerance_unit')!='Da':raise ValueError('Relative windows require explicit Da uncertainty')
        basis=observation.get('mass_basis','monoisotopic')
        if basis not in {'monoisotopic','average'} or kind=='delta_mono' and basis!='monoisotopic':raise ValueError('Relative mass basis disagrees with its kind')
        return [{'kind':'delta','value':value,'tolerance':tolerance,'basis':basis,'label':'explicit_relative_window'}]
    if kind=='peak_pair':
        check_keys(observation,{'kind','peaks','error_model','reference','reference_chain_ids'},'peak pair')
        peaks=observation.get('peaks')
        if not isinstance(peaks,list) or len(peaks)!=2:raise ValueError('A peak pair requires exactly two measured peaks')
        options=[]
        for first in hypotheses(peaks[0]):
            for second in hypotheses(peaks[1]):
                pair=peak_pair(first,second,observation.get('error_model','bounded'))
                options.append({'kind':'pair','pair':pair,'first':first,'second':second,'basis':pair['mass_basis'],
                                'value':pair['delta_element_mass'],'tolerance':pair['tolerance_da'],
                                'label':first.assumption+'/'+second.assumption})
        return options
    return [{'kind':'absolute','coordinate':x,'basis':x.mass_basis,'label':x.assumption} for x in hypotheses(observation)]

def event_summary(event):
    return {k:deepcopy(v) for k,v in event.items() if k not in {'product_state','contract','input_binding_signature'}}

def diagnose_state_case(case):
    check_keys(case,{'case_id','reaction','context','intended_product','starting_material','observations','max_events','max_nodes','rule_ids','input_text'},'state case')
    if not isinstance(case.get('observations'),list) or not case['observations']:raise ValueError('Provide measured observations')
    bundle,version=load_runtime()
    target=normalize_molecule(case['intended_product'])
    starting=normalize_molecule(case['starting_material']) if case.get('starting_material') else None
    context=prepare_context(case,target,starting)
    depth=case.get('max_events',2);budget=case.get('max_nodes',100000)
    search=enumerate_states(bundle,target,context,depth,budget,case.get('rule_ids'))
    roots=[('intended_product',search)]
    if starting is not None:
        roots.append(('starting_material',enumerate_states(bundle,starting,context,depth,budget,case.get('rule_ids'))))
    observations=[];component_cache={}
    for observation in case['observations']:
        options=observation_options(observation);refs=reference_components(case,observation,target,starting)
        matches={};conflicts={};path_signatures={}
        for origin,run in roots:
            for node in run['nodes']:
                state_key=fingerprint(node['state'])
                if state_key not in component_cache:component_cache[state_key]=describe(node['state'])
                for component in component_cache[state_key]:
                    q=component['intrinsic_charge'];paths=None
                    for option in options:
                        predicted_mass=formula_mass(component['formula'],option['basis'])
                        for ref in ([None] if option['kind']=='absolute' else refs):
                            if option['kind']=='absolute':
                                coordinate=option['coordinate']
                                if q!=coordinate.intrinsic_charge:continue
                                result=coordinate.match(predicted_mass)
                            else:
                                if option['kind']=='delta' and (q or ref['intrinsic_charge']):continue
                                if option['kind']=='pair' and (q!=option['second'].intrinsic_charge or ref['intrinsic_charge']!=option['first'].intrinsic_charge):continue
                                predicted=predicted_mass-formula_mass(ref['formula'],option['basis'])
                                error=option['value']-predicted
                                result={'match':abs(error)<=option['tolerance'],'predicted':predicted,'error_da':error,'error_ppm':None}
                            if not result['match']:continue
                            if paths is None:paths=[event_summary(e) for e in node['events']]
                            blocked=any(e['condition']=='conflict' for e in paths)
                            bucket=conflicts if blocked else matches
                            key=(physical_fingerprint(component['state']),option['label'],tuple(ref['chain_ids']) if ref else ())
                            if key not in bucket:
                                bucket[key]={'component':component,'predicted':result['predicted'],'error_da':result['error_da'],
                                    'error_ppm':result.get('error_ppm'),'mass_basis':option['basis'],'ion_assumption':option['label'],
                                    'reference_chain_ids':ref['chain_ids'] if ref else None,'origins':[],'paths':[],
                                    'condition':'unknown' if any(e['condition']=='unknown' for e in paths) else 'compatible',
                                    'evidence':'mass_match_only'}
                                if option['kind']=='pair':bucket[key]['propagated_error']=option['pair']
                            hit=bucket[key]
                            if origin not in hit['origins']:hit['origins'].append(origin)
                            signatures=path_signatures.setdefault((blocked,key),set())
                            signature=encoded([(e['branch_id'],e['sites'],e.get('reacting_endpoint')) for e in paths])
                            if signature not in signatures:hit['paths'].append(paths);signatures.add(signature)
                            if all(e['condition']=='compatible' for e in paths):hit['condition']='compatible'
        accepted=list(matches.values());excluded=list(conflicts.values())
        accepted.sort(key=lambda x:(min(len(p) for p in x['paths']),abs(x['error_da']),encoded(x['component']['formula'])))
        observations.append({'input':observation,'candidates':accepted,'condition_conflicts':excluded,
                             'candidate_count':len(accepted),'conflict_count':len(excluded),
                             'status':'matched_hypotheses' if accepted else 'no_match_in_searched_pool',
                             'ion_hypothesis_count':len(options),
                             'source_leads':mass_source_leads(bundle,observation,options,refs,context,target)})
    scopes=[{'origin':origin,**run['scope']} for origin,run in roots]
    return {'engine':'codex-state-core-1','library_hash':version,'case':case,'reference_components':describe(target),
            'starting_components':describe(starting) if starting else [],'observations':observations,'search_scopes':scopes,
            'complete_for_compiled_transition_pool':all(x['complete_for_compiled_transition_pool'] for x in scopes),
            'full_book_capabilities_complete':all(x['full_book_capabilities_complete'] for x in scopes),
            'complete_for_practical_mass_scope':all(not x['truncated'] and not x['unsupported_transitions'] for x in scopes),
            'scope_grade':'practical_mass_families_and_supported_state_transitions',
            'real_sample':'not_verified'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path,nargs='?');parser.add_argument('--text');parser.add_argument('--output',type=Path);parser.add_argument('--report',type=Path)
    parser.add_argument('--max-events',type=int);parser.add_argument('--max-nodes',type=int);parser.add_argument('--lookup');parser.add_argument('--stats',action='store_true')
    parser.add_argument('--find')
    args=parser.parse_args()
    try:
        if args.lookup or args.stats or args.find:
            from knowledge_index import lookup,find_records
            bundle,version=load_runtime()
            result={'library_hash':version,'counts':{'rules':len(bundle['rules']),'branches':sum(len(r['branches']) for r in bundle['rules']),'appendix':len(bundle['appendix_dispositions'])}} if args.stats else {'query':args.find,'records':find_records(args.find)} if args.find else lookup(args.lookup)
            if result is None:result=next((r for r in bundle['appendix_dispositions'] if r['rule_id']==args.lookup),None)
            if result is None:raise ValueError('Unknown rule/appendix ID')
        elif bool(args.case)==bool(args.text):raise ValueError('Provide one JSON case or one --text input')
        elif args.text:
            from natural_case import parse_text
            parsed=parse_text(args.text)
            if parsed['status']=='ready':
                case=parsed['case']
                if args.max_events is not None:case['max_events']=args.max_events
                if args.max_nodes is not None:case['max_nodes']=args.max_nodes
                result=diagnose_state_case(case)
            else:
                observation=parsed.get('parsed',{}).get('observation')
                if observation and observation.get('kind') in {'delta_mono','delta_approx'}:
                    from source_leads import lookup_delta
                    from case_input import read_sequence
                    sequence=parsed['parsed'].get('sequence')
                    conditions=parsed['parsed'].get('conditions',{})
                    result=lookup_delta(load_runtime()[0],observation['value'],observation.get('tolerance') if observation.get('tolerance_unit')=='Da' else None,
                                        read_sequence(sequence) if sequence else None,conditions.get('stage_tags',[]),conditions.get('stage_tags_absent',[]))
                    result.update(missing=parsed['missing'],input_text=args.text)
                else:result=parsed
        else:
            case=json.loads(args.case.read_text(encoding='utf-8-sig'))
            if args.max_events is not None:case['max_events']=args.max_events
            if args.max_nodes is not None:case['max_nodes']=args.max_nodes
            result=diagnose_state_case(case)
    except (ValueError,KeyError,TypeError) as error:parser.error(str(error))
    text=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)
    if args.output:args.output.write_text(text+'\n',encoding='utf-8')
    else:print(text)
    if args.report:
        from state_report import render_state_report
        if result.get('engine'):report=render_state_report(result)
        elif result.get('lookup_mode'):
            report='# 质量家族粗查\n\n'+result['precision_note']+'。\n\n'
            report+='\n'.join('- '+x['id']+'：'+str(x['name'])+'；'+x['context_status'] for x in result['source_leads'])+'\n\n'
            report+='完整状态计算尚需：'+'、'.join(result.get('missing',[]))+'。\n'
        elif result.get('missing'):report='# 尚需补充的信息\n\n'+'\n'.join('- '+x for x in result['missing'])+'\n'
        else:report='# 知识库查询\n\n```json\n'+text+'\n```\n'
        args.report.write_text(report,encoding='utf-8')

if __name__=='__main__':main()
