"""Formal Codex/DSH fused CLI using one canonical library and one structural executor."""
import argparse,json
from pathlib import Path
from fusion_bridge import diagnose_fused_case,adapt_dsh_case
from knowledge_source import load_runtime
from natural_case import parse_text
from source_leads import lookup_delta
from state_report import render_state_report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',nargs='?',type=Path)
    p.add_argument('--text');p.add_argument('--output',type=Path);p.add_argument('--report',type=Path)
    p.add_argument('--max-events',type=int);p.add_argument('--max-nodes',type=int)
    p.add_argument('--lookup');p.add_argument('--find');p.add_argument('--stats',action='store_true')
    a=p.parse_args()
    try:
        if a.stats or a.lookup or a.find:
            from knowledge_index import lookup,find_records
            bundle,version=load_runtime()
            result={'engine':'codex-dsh-fused-1','library_hash':version,'counts':{'rules':len(bundle['rules']),
                'branches':sum(len(r['branches']) for r in bundle['rules']),'appendix':len(bundle['appendix_dispositions'])}} if a.stats else {'records':find_records(a.find)} if a.find else lookup(a.lookup)
            if result is None:result=next((r for r in bundle['appendix_dispositions'] if r['rule_id']==a.lookup),None)
            if result is None:raise ValueError('Unknown source ID')
        else:
            if bool(a.case)==bool(a.text):raise ValueError('Provide one JSON case or --text')
            if a.text:
                parsed=parse_text(a.text)
                if parsed['status']!='ready':
                    obs=parsed['parsed'].get('observation')
                    if obs and obs['kind'] in {'delta_mono','delta_approx'}:
                        from case_input import read_sequence
                        c=parsed['parsed']['conditions'];seq=parsed['parsed'].get('sequence')
                        result=lookup_delta(load_runtime()[0],obs['value'],obs.get('tolerance') if obs.get('tolerance_unit')=='Da' else None,
                            read_sequence(seq) if seq else None,c['stage_tags'],c['stage_tags_absent'])
                        result.update(missing=parsed['missing'],input_text=a.text)
                    else:result=parsed
                    case=None
                else:case=parsed['case']
            else:
                case=json.loads(a.case.read_text(encoding='utf-8-sig'))
                if 'sequence' in case:case=adapt_dsh_case(case)
            if case is not None:
                if a.max_events is not None:case['max_events']=a.max_events
                if a.max_nodes is not None:case['max_nodes']=a.max_nodes
                result=diagnose_fused_case(case)
    except (ValueError,KeyError,TypeError) as error:p.error(str(error))
    output=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
    if a.output:a.output.write_text(output,encoding='utf-8')
    else:print(output,end='')
    if a.report:
        if 'search_scopes' in result:report=render_state_report(result)
        else:report='# 多肽质量来源粗查\n\n'+output+'\n未进行完整结构搜索；缺项：'+str(result.get('missing',[]))
        a.report.write_text(report,encoding='utf-8')
if __name__=='__main__':main()
