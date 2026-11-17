"""Independent small-system endpoint enumeration and tolerance/irrelevant-context invariants."""
import itertools,json
from pathlib import Path
from pyteomics import mass
from state_diagnosis import diagnose_state_case
ROOT=Path(__file__).resolve().parents[1]
effects={'ox':('CH09-R001',mass.calculate_mass(formula='O')),
    'delete':('CH01-R002',-mass.calculate_mass(sequence='S')+mass.calculate_mass(formula='H2O')),
    'nform':('CH07-R013',mass.calculate_mass(formula='CO')),
    'sform':('CH07-R013',mass.calculate_mass(formula='CO'))}
expected={}
for depth in range(3):
    for actions in itertools.permutations(effects,depth):
        if 'delete' in actions and 'sform' in actions and actions.index('delete')<actions.index('sform'):continue
        delta=sum(effects[a][1] for a in actions)
        if 'delete' in actions and 'sform' in actions:delta-=effects['sform'][1]
        expected.setdefault(round(delta,6),set()).add(tuple(sorted(effects[a][0] for a in actions)))
observations=[{'kind':'delta_mono','value':d,'tolerance':0.00001,'tolerance_unit':'Da'} for d in expected]
case={'case_id':'independent-sequence-combinations','intended_product':{'sequence':'ACSA','nterm':'free','cterm':'acid'},
      'reaction':{'description':'independent endpoint enumeration'},'observations':observations,'max_events':2,'rule_ids':sorted({x[0] for x in effects.values()})}
result=diagnose_state_case(case);checks=[]
for delta,o in zip(expected,result['observations'],strict=True):
    actual={tuple(sorted(e['rule_id'] for e in path)) for h in o['candidates'] for path in h['paths']}
    checks.append({'expected_rules':sorted(expected[delta]),'actual_rules':sorted(actual),'pass':actual==expected[delta],'delta':delta})
wide=json.loads(json.dumps(case));wide['reaction']['stage_tags']=['unrelated_fact'];
for o in wide['observations']:o['tolerance']=0.001
expanded=diagnose_state_case(wide)
subset=True
for old,new in zip(result['observations'],expanded['observations'],strict=True):
    old_ids={tuple(sorted(e['branch_id'] for e in p)) for h in old['candidates'] for p in h['paths']}
    new_ids={tuple(sorted(e['branch_id'] for e in p)) for h in new['candidates'] for p in h['paths']}
    subset&=old_ids<=new_ids
report={'expected_endpoints':len(expected),'checks':checks,'tolerance_and_irrelevant_context_monotonicity':subset,
        'pass':all(c['pass'] for c in checks) and subset,'oracle':'independent Cys oxidation/Ser deletion/N-formyl/Ser O-formyl ordered choices; deleting Ser removes its O-formyl; Pyteomics formulas, no transform reused'}
(ROOT/'rebuild/practical_independent_oracles.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False))
