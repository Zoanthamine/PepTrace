"""Five public exercises, scored from the fused CLI output rather than reference prose."""
import json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'validation/five_results'
def run(case,path,max_events=None):
    args=[sys.executable,'-B',str(ROOT/'scripts/diagnose.py'),str(case),'--output',str(path)]
    if max_events:args+=['--max-events',str(max_events)]
    subprocess.run(args,cwd=ROOT,check=True);return json.loads(path.read_text(encoding='utf-8'))
def score(i,result):
    o=result['observations'][0];paths=[p for h in o['candidates'] for p in h['paths']]
    events=[e for p in paths for e in p]
    checks={'mass_match':bool(o['candidates']),'not_truncated':not any(s['truncated'] for s in result['search_scopes']),
        'independent_formula':result['fusion_audit']['whole_product_formula_checks_pass']}
    if i==1:
        checks.update(n_terminal_phe=any(e['branch_id']=='CH10-R003:001' for e in events),
            lys_sidechain_phe=any(e['branch_id']=='CH10-R003:002' for e in events))
    elif i==2:
        mods=[m for h in o['candidates'] for c in h['component']['state']['chains'] for m in c['peptide']['modifications']]
        for site in (2,6,'C'):checks['methyl_ester_'+str(site)]=any(m['site']==site and m.get('delta_formula')=={'C':1,'H':2} and '甲酯' in m.get('name','') for m in mods)
    elif i in (3,4):
        checks['trisulfide']=any(b['kind']=='trisulfide' for h in o['candidates'] for b in h['component']['state']['bonds'])
        if i==3:checks['DCC_path']=any(e['branch_id']=='CH13-R014:000' for e in events)
        else:
            checks['DCC_excluded']=any(e['branch_id']=='CH13-R014:000' and 'dcc' in e['conflicting_conditions'] for h in o['condition_conflicts'] for p in h['paths'] for e in p)
            checks['DCC_not_retained']=not any(e['branch_id']=='CH13-R014:000' for e in events)
    else:
        for site in (2,3,4):checks['tert_butyl_'+str(site)]=any(e['rule_id']=='CH14-R020' and f'chain1:{site}' in e['sites'] for e in events)
        checks['acid_not_guessed']=any('strong_acid_workup' in e['missing_conditions'] for e in events if e['rule_id']=='CH14-R020')
    return checks
def main():
    OUT.mkdir(parents=True,exist_ok=True);rows=[]
    for i in range(1,6):
        start=time.perf_counter();r=run(ROOT/f'examples/five_cases/case-{i}.json',OUT/f'result-{i}.json')
        checks=score(i,r);o=r['observations'][0]
        row={'case':i,'checks':checks,'pass':all(checks.values()),'seconds':time.perf_counter()-start,
            'candidates':o['candidate_count'],'conflicts':o['conflict_count'],'nodes':sum(s['nodes'] for s in r['search_scopes'])}
        rows.append(row);print(json.dumps(row,ensure_ascii=False),flush=True)
    near=json.loads((ROOT/'examples/five_cases/case-3.json').read_text(encoding='utf-8'));near['max_events']=1
    for p in near['observations'][0]['peaks']:p['tolerance']=.005
    p=OUT/'resolution-input.json';p.write_text(json.dumps(near),encoding='utf-8');n=run(p,OUT/'resolution-result.json')
    ncheck=any(b['kind']=='thiosulfonate' for h in n['observations'][0]['candidates'] for b in h['component']['state']['bonds'])
    charged=json.loads((ROOT/'examples/five_cases/case-5.json').read_text(encoding='utf-8'));charged['max_events']=1
    charged['observations']=[{'kind':'mz','value':785.35,'tolerance':.1,'tolerance_unit':'Da','ion_hypotheses':[
        {'charge':1,'ions':{'H':1}},{'charge':1,'ions':{},'intrinsic_charge':1}]}]
    p=OUT/'charge-input.json';p.write_text(json.dumps(charged),encoding='utf-8');c=run(p,OUT/'charge-result.json')
    ccheck=any(h['component']['intrinsic_charge']==1 and any(e['branch_id']=='CH14-R020:001' for p in h['paths'] for e in p) for h in c['observations'][0]['candidates'])
    report={'type':'public exercises; not blinded or real-sample accuracy','rows':rows,
        'resolution_sensitivity':ncheck,'fixed_charge_sensitivity':ccheck,'pass':all(x['pass'] for x in rows) and ncheck and ccheck}
    (ROOT/'validation/five_cases.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if not report['pass']:raise SystemExit(1)
if __name__=='__main__':main()
