"""Cold/warm baseline first, fixed practical thresholds second; preserve full search counts."""
import json,subprocess,sys,time,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sample(seq):return {'case_id':'performance-'+seq,'intended_product':{'sequence':seq,'nterm':'free','cterm':'acid'},
    'reaction':{'description':'performance baseline','stage_tags':['oxidation']},'max_events':2,'max_nodes':100000,
    'observations':[{'kind':'delta_mono','value':15.99491461957,'tolerance':0.001,'tolerance_unit':'Da'}]}
def main():
    from state_diagnosis import diagnose_state_case
    rows=[]
    for seq in ('M','GC','MELK'):
        case=sample(seq);warm=[];counts=[]
        for _ in range(3):
            before=time.perf_counter();result=diagnose_state_case(case);warm.append(time.perf_counter()-before)
            scopes=result['search_scopes'];counts.append({'nodes':sum(s['nodes'] for s in scopes),'candidates':result['observations'][0]['candidate_count'],
                                                         'truncated':any(s['truncated'] for s in scopes)})
        file=ROOT/'rebuild'/('benchmark-'+seq+'.json');file.write_text(json.dumps(case),encoding='utf-8')
        cold=[]
        for _ in range(2):
            before=time.perf_counter();run=subprocess.run([sys.executable,'-B','tools/state_diagnosis.py',str(file),'--output',str(ROOT/'rebuild/benchmark-result.json')],
                                                         cwd=ROOT,capture_output=True,text=True,encoding='utf-8',timeout=60)
            if run.returncode:raise RuntimeError(run.stderr)
            cold.append(time.perf_counter()-before)
        rows.append({'sequence':seq,'length':len(seq),'event_budget':2,'node_budget':100000,'warm_seconds':warm,'cold_seconds':cold,'counts':counts})
    # Freeze a threshold after observing the baseline, then perform a fresh verification run.
    threshold={'warm_seconds':10.0,'cold_seconds':10.0}
    verification=[]
    for seq in ('M','GC','MELK'):
        before=time.perf_counter();result=diagnose_state_case(sample(seq));elapsed=time.perf_counter()-before
        verification.append({'sequence':seq,'seconds':elapsed,'pass':elapsed<=threshold['warm_seconds'] and not any(s['truncated'] for s in result['search_scopes'])})
    data={'baseline':rows,'fixed_after_baseline_threshold':threshold,'fresh_verification':verification,
          'pass':all(x['pass'] for x in verification),'scope':'short practical sequences, all compiled rules, two applications; not a hardware stress test'}
    (ROOT/'rebuild/practical_benchmark.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'pass':data['pass'],'threshold':threshold,'workloads':[{'seq':r['sequence'],'warm_median':statistics.median(r['warm_seconds']),
                                                                          'cold_median':statistics.median(r['cold_seconds']),'counts':r['counts'][0]} for r in rows]},ensure_ascii=False))
if __name__=='__main__':main()
