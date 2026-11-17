"""Verify package integrity and runnable behavior; --full executes the meaningful regression gates."""
import argparse,importlib.metadata,json,re,subprocess,sys,time
from pathlib import Path
from check_release import check
ROOT=Path(__file__).resolve().parents[1];ENGINE=ROOT/'engine'
def run(command,cwd,label):
    start=time.perf_counter();p=subprocess.run([sys.executable,'-B',*command],cwd=cwd,capture_output=True,text=True,encoding='utf-8')
    log=ROOT/'validation'/f'{label}.log';log.write_text(p.stdout+p.stderr,encoding='utf-8')
    if p.returncode:raise RuntimeError(f'{label} failed; see validation/{label}.log')
    return {'label':label,'seconds':time.perf_counter()-start,'returncode':p.returncode,'output':p.stdout+p.stderr}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--benchmark',action='store_true');a=parser.parse_args()
    (ROOT/'validation').mkdir(exist_ok=True);integrity=check();rows=[]
    smoke=run([str(ROOT/'scripts/diagnose.py'),'--stats'],ROOT,'portable-stats')
    stats=json.loads(smoke['output']);assert stats['engine']=='codex-dsh-fused-1' and stats['counts']['branches']==308
    run([str(ROOT/'scripts/diagnose.py'),'--text','序列 MELK，N端游离，C端羧酸，比理论多15.994915 Da，容差0.001 Da，用DMSO处理过','--max-events','1','--output',str(ROOT/'validation/smoke.json')],ROOT,'portable-nl')
    r=json.loads((ROOT/'validation/smoke.json').read_text(encoding='utf-8'))
    assert any(e['rule_id']=='CH09-R029' for h in r['observations'][0]['candidates'] for p in h['paths'] for e in p)
    core_count=0;package_count=0
    if a.full:
        t=run(['-m','unittest','discover','-s','tools','-p','test_*.py'],ENGINE,'regression')
        core_count=int(re.search(r'Ran (\d+) tests?',t['output'])[1]);rows.append({'label':'core_regression','tests':core_count,'pass':True})
        t=run(['-m','unittest','discover','-s','tests'],ROOT,'package-regression')
        package_count=int(re.search(r'Ran (\d+) tests?',t['output'])[1]);rows.append({'label':'package_regression','tests':package_count,'pass':True})
        for script,resultpath in [('branch_acceptance.py','rebuild/branch_acceptance.json'),
            ('practical_independent_oracles.py','rebuild/practical_independent_oracles.json'),('run_practical_holdouts.py','rebuild/holdouts/report.json')]:
            run([str(ENGINE/'tools'/script)],ENGINE,Path(script).stem)
            result=json.loads((ENGINE/resultpath).read_text(encoding='utf-8'))
            passed=result['failed']==0 if 'failed' in result else result['pass'];assert passed
            rows.append({'label':Path(script).stem,'pass':passed,'total':result.get('total')})
        run([str(ROOT/'scripts/check_five_cases.py')],ROOT,'five-cases')
        five=json.loads((ROOT/'validation/five_cases.json').read_text(encoding='utf-8'));assert five['pass'];rows.append({'label':'five_public_cases','pass':True,'total':5})
    if a.benchmark:
        run([str(ENGINE/'tools/practical_benchmark.py')],ENGINE,'benchmark')
        b=json.loads((ENGINE/'rebuild/practical_benchmark.json').read_text(encoding='utf-8'));assert b['pass'];rows.append({'label':'benchmark','pass':True})
    versions={}
    for name in ['pyteomics','jsonschema','rdkit','numpy','pyopenms']:
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]='not_installed'
    report={'version':'1.0.0','integrity':integrity,'library_hash':stats['library_hash'],'full':a.full,
        'core_regression_tests':core_count,'package_regression_tests':package_count,'checks':rows,
        'environment':{'python':sys.version.split()[0],'packages':versions},'pass':True,'real_sample':'not_verified'}
    (ROOT/'validation/acceptance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
