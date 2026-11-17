"""Validate canonical sources and build all compatibility exports plus a frozen runtime."""
from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
from pyteomics import mass
from jsonschema import Draft202012Validator
from rule_contracts import compile_contract, validate_contract, legacy_route
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--source-root',type=Path,help='Optional directory containing doc_XX/source.pdf for original-document hash rechecks')
args=parser.parse_args()
source_checks={}
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def content_hash(v):return hashlib.sha256(json.dumps(v,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def check_formula(f,absolute=False):
    assert isinstance(f,dict)
    assert all(e in constants['monoisotopic_atoms'] and type(n) is int and n!=0 for e,n in f.items()),f
    assert not absolute or all(n>0 for n in f.values()),f
def fm(f):return sum(constants['monoisotopic_atoms'][e]*n for e,n in f.items())

constants=read(ROOT/'knowledge/mass_constants.json');components=read(ROOT/'knowledge/components.json')
rules=read(ROOT/'knowledge/rules.json');records=rules['records']
schema=read(ROOT/'knowledge/schema.json')
for definition,value in [('components',components),('rules',rules)]:
    Draft202012Validator({'$ref':'#/$defs/'+definition,'$defs':schema['$defs']}).validate(value)
ids=[r['rule_id'] for r in records];assert len(ids)==len(set(ids))
assert sum(r['legacy_record'] is not None for r in records)==rules['migration']['rule_source_count']
all_branches=[b for r in records for b in r['branches']]
assert len({b['branch_id'] for b in all_branches})==len(all_branches)
for row in records:
    row['execution_route']=legacy_route(row)
    for branch in row['branches']:
        recipe=branch['recipe'];selector=recipe.get('site_selector',recipe.get('bond_selector',{}));transform=recipe['transform']
        for field in ('required_conditions','required_any_conditions','required_features','supporting_any_conditions'):
            values=recipe.get(field,[])
            assert isinstance(values,list) and all(isinstance(x,str) and x for x in values),(branch['branch_id'],field)
        group_refs=set(selector.get('groups',[]))
        if selector.get('group'):group_refs.add(selector['group'])
        for field in ('product_group','new_group'):
            if transform.get(field):group_refs.add(transform[field])
        assert group_refs<=set(components['covalent_groups']),(branch['branch_id'],group_refs-set(components['covalent_groups']))
        residue_refs=set()
        for field in ('residues','requires_residues','exclude_residues','next_residues','motif'):residue_refs.update(selector.get(field,[]))
        if transform.get('residue'):residue_refs.add(transform['residue'])
        if transform.get('residues_by_input'):residue_refs.update(transform['residues_by_input']);residue_refs.update(transform['residues_by_input'].values())
        assert residue_refs<=set(components['residues']),(branch['branch_id'],residue_refs-set(components['residues']))
        for field in ('delta_formula','replacement_formula','cap_formula','bond_delta','migrated_formula'):
            if field in transform:check_formula(transform[field])
        branch['contract']=compile_contract(branch['recipe'])
        validate_contract(branch['contract'])
checks=[]
for name,r in components['residues'].items():
    a,b=r['residue_formula'],r['free_amino_acid_formula'];check_formula(a,True);check_formula(b,True)
    assert {e:b.get(e,0)-a.get(e,0) for e in set(a)|set(b) if b.get(e,0)!=a.get(e,0)}=={'H':2,'O':1}
    checks.extend([{'kind':'residue','name':name,'expected':mass.calculate_mass(composition=a),'actual':fm(a)},
                   {'kind':'free_amino_acid','name':name,'expected':mass.calculate_mass(composition=b),'actual':fm(b)}])
for category,key in [('covalent_groups','net_formula'),('neutral_molecules','formula')]:
    for name,r in components[category].items():
        f=r[key];check_formula(f,category=='neutral_molecules')
        checks.append({'kind':category,'name':name,'expected':mass.calculate_mass(composition=f),'actual':fm(f)})
for c in checks:assert abs(c['expected']-c['actual'])<1e-5,c
for row in records:
    t=row['template']
    overrides=[x['content'] for x in row['review_layers'] if x['kind']=='source_override']
    row['effective_source']=overrides[-1] if overrides else (row['legacy_record'] or {}).get('source')
    if not t:continue
    assert t['rule_id']==row['rule_id'] and t['required_conditions'] is not None
    assert t['reference_state'] and t['evidence_level'] and t['source_refs']
    if 'delta_formula' in t['transform']:check_formula(t['transform']['delta_formula'])
    source=rules['pack_metadata'][row['compatibility_pack']]['source']
    if source.get('doc_id'):
        if args.source_root is not None:
            p=args.source_root/source['doc_id']/'source.pdf'
            assert source['sha256']==hashlib.sha256(p.read_bytes()).hexdigest(),row['rule_id']
            source_checks[source['doc_id']]='hash_verified'
        else:source_checks[source['doc_id']]='original_not_bundled; canonical facts rebuilt without a new visual/source review'
    for key in ('group','new_group','required_group'):
        val=t['transform'].get(key)
        if val:assert val in components['covalent_groups'],(row['rule_id'],val)
for row in rules['appendix_dispositions']:
    assert set(row['legacy_rule_ids'])<=set(ids),row['rule_id']
    assert set(row['new_rule_ids'])<=set(ids),row['rule_id']
for row in records:
    d=row['disposition']
    if d:assert set(d['new_rule_ids'])<=set(ids)

# Rebuild every compatibility mass table from the canonical formula, never copy stored floats.
compat={'schema_version':2,'generated_from':'knowledge/components.json; do not edit',
        'residues':{},'covalent_deltas':{},'neutral_molecules':{},'sequences_free_ends':{}}
for name,r in components['residues'].items():
    compat['residues'][name]={'residue':{'formula':r['residue_formula'],'monoisotopic_da':fm(r['residue_formula'])},
                            'free_amino_acid':{'formula':r['free_amino_acid_formula'],'monoisotopic_da':fm(r['free_amino_acid_formula'])}}
for name,r in components['covalent_groups'].items():
    compat['covalent_deltas'][name]={'formula':r['net_formula'],'reference':r['reference_state'],'monoisotopic_da':fm(r['net_formula'])}
for name,r in components['neutral_molecules'].items():
    compat['neutral_molecules'][name]={'formula':r['formula'],'monoisotopic_da':fm(r['formula'])}
for seq,label in components['sequence_examples'].items():
    f={'H':2,'O':1}
    for name in seq:
        for e,n in components['residues'][name]['residue_formula'].items():f[e]=f.get(e,0)+n
    compat['sequences_free_ends'][seq]={'formula':f,'monoisotopic_da':fm(f),'description':label}
write(ROOT/'多肽质量组件.json',compat)
write(ROOT/'规则库_总.json',{**rules['legacy_book_metadata'],'generated_from':'knowledge/rules.json; do not edit',
                             'rules':[r['legacy_record'] for r in records if r['origin']=='book']})
write(ROOT/'SUPPLEMENTS.json',{**rules['legacy_supplement_metadata'],'generated_from':'knowledge/rules.json; do not edit',
                              'items':[r['legacy_record'] for r in records if r['origin']=='supplement']})
for name,meta in rules['pack_metadata'].items():
    write(ROOT/'rules'/name,{**meta,'generated_from':'knowledge/rules.json; do not edit',
                            'rules':[r['template'] for r in records if r['compatibility_pack']==name]})
appendix=read(ROOT/'knowledge/appendix.json')
adapters=read(ROOT/'knowledge/execution_adapters.json')
assert set(adapters['relative_profiles'])<=set(ids)
assert set(adapters['site_rules'])<=set(ids)
write(ROOT/'appendix_I_recovered.json',appendix['source_rows'])
write(ROOT/'appendix_I_element_account.json',appendix['element_accounts'])
runtime={'schema_version':'codex-runtime-1','components':components,'mass_constants':constants,
         'rules':records,'appendix_dispositions':rules['appendix_dispositions'],'appendix':appendix,
         'pack_metadata':rules['pack_metadata'],'execution_adapters':adapters}
runtime['continuation_paths']=read(ROOT/'knowledge/continuation_paths.json')
runtime['release_scope']=read(ROOT/'knowledge/release_scope.json')
assert {p['next_rule_id'] for p in runtime['continuation_paths']['paths']}<=set(ids)
write(ROOT/'runtime/knowledge.json',runtime)
sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'knowledge').glob('*.json'))}
exports=[ROOT/'规则库_总.json',ROOT/'SUPPLEMENTS.json',ROOT/'多肽质量组件.json',
         ROOT/'appendix_I_recovered.json',ROOT/'appendix_I_element_account.json',*(ROOT/'rules').glob('*.json')]
manifest={'schema_version':'codex-build-manifest-1','canonical_root':'knowledge',
          'source_sha256':sources,'logical_sha256':content_hash(runtime),
          'runtime_sha256':hashlib.sha256((ROOT/'runtime/knowledge.json').read_bytes()).hexdigest(),
          'compatibility_exports':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(exports)},
          'compiler_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'counts':{'rules':len(records),'templates':sum(bool(r['template']) for r in records),'branches':len(all_branches),
                    'appendix':len(rules['appendix_dispositions']),'mass_components_checked':len(checks)},
          'state_constraints':'compiled practical state contracts; exhaustive atom/geometry chemistry outside release scope',
          'branch_acceptance':'see validation/acceptance.json', 'original_source_recheck':source_checks,'real_sample':'not_verified'}
write(ROOT/'runtime/manifest.json',manifest)
write(ROOT/'runtime/legacy_id_map.json',{r['rule_id']:r['execution_route'] for r in records})
write(ROOT/'rebuild/r02_contract_validation.json',{'schema':'state-rule-contract-1','contracts':len(all_branches),
    'families':dict(Counter(b['contract']['transform_family'] for b in all_branches)),
    'legacy_routes':len(records),'all_contracts_valid':True,
    'boundary':'Compiled contracts describe existing recipes; runtime state execution and all source-model gaps remain separate gates.'})
write(ROOT/'rebuild/r02_mass_validation.json',{'checks':checks,'all_pass':True,'independent_backend':'pyteomics.mass'})
print(json.dumps(manifest['counts'],ensure_ascii=False))
