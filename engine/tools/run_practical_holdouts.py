"""Solve frozen questions in separate processes, then reveal the independent answer key."""
import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'rebuild/holdouts'
registration=json.loads((OUT/'registration.json').read_text(encoding='utf-8'));results=[]
for i,entry in enumerate(registration['questions'],1):
    p=OUT/entry['file'];assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
    run=subprocess.run([sys.executable,'-B','tools/state_diagnosis.py',str(p),'--output',str(OUT/f'result-{i}.json')],
                       cwd=ROOT,capture_output=True,text=True,encoding='utf-8',timeout=60)
    if run.returncode:raise RuntimeError(run.stderr)
    results.append(json.loads((OUT/f'result-{i}.json').read_text(encoding='utf-8')))
# The solver processes above never receive the answers path or file contents.
answers=json.loads((OUT/'answers.json').read_text(encoding='utf-8'));checks=[]
for answer,result in zip(answers,results,strict=True):
    found={x['id'] for o in result['observations'] for x in o['source_leads']}
    checks.append({'case_id':answer['case_id'],'expected_appendix_id':answer['expected_appendix_id'],
                   'pass':answer['expected_appendix_id'] in found,'candidate_counts':[o['candidate_count'] for o in result['observations']],
                   'truncated':any(s['truncated'] for s in result['search_scopes'])})
report={'seed':registration['seed'],'pool_size':len(registration['pool']),'tests':checks,'pass':all(x['pass'] and not x['truncated'] for x in checks),
        'scope':'three random literature mass-family retrieval holdouts; not full-library or real-sample structural identification accuracy',
        'answer_isolation':'solver received frozen knowledge and questions only; answers read after all solver processes exited'}
(OUT/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False))
