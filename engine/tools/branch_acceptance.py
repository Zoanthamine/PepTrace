"""Practical per-branch execution/condition checks with an independent Pyteomics mass oracle."""
from copy import deepcopy
import json
from pathlib import Path
from pyteomics import mass
from knowledge_source import load_runtime
from structure_state import normalize_molecule
from state_transformations import linear_transitions,bridge_transitions
from carbonyl_components import normalize_carbonyl_donors

ROOT=Path(__file__).resolve().parents[1]
def merge(*fs):
    out={}
    for f in fs:
        for e,n in f.items():out[e]=out.get(e,0)+n
    return {e:n for e,n in out.items() if n}
def negate(f):return {e:-n for e,n in f.items()}
def m(f):return mass.calculate_mass(composition=f)

def fixture(bundle,branch):
    r=branch['recipe'];s=r.get('site_selector',{});t=r['transform'];kind=t['kind'];comp=bundle['components']
    value=next((x for x in s.get('residues',[]) if x not in s.get('exclude_residues',[])),None)
    if value is None:value=next(x for x in ['A','L','G','S','C'] if x not in s.get('exclude_residues',[]))
    motif=s.get('motif') or [value];position=s.get('position','any')
    seq=list(motif)+['A','G'] if position=='start' else ['A']+list(motif) if position=='end' else ['A']+list(motif)+['G']
    index=0 if position=='start' else 1
    if s.get('previous_required') and index==0:seq.insert(0,'A');index+=1
    if s.get('next_residues'):
        if index+1==len(seq):seq.append(s['next_residues'][0])
        else:seq[index+1]=s['next_residues'][0]
    if s.get('next_is_end'):seq=seq[:index+2]
    if t.get('previous_residue'):
        if index==0:seq.insert(0,t['previous_residue']);index+=1
        else:seq[index-1]=t['previous_residue']
    for residue in s.get('requires_residues',[]):
        if residue not in seq:seq.append(residue)
    if s.get('sequence_length'):
        seq=[value]+['A']*(s['sequence_length']-1);index=0
    peptide={'sequence':seq,'nterm':s.get('nterm','free'),'cterm':s.get('cterm','acid'),'modifications':[]}
    group=s.get('group') or (s.get('groups') or [None])[0]
    if group:
        site=s.get('group_site',index+1)
        if group in {'NCA','NChloroformylAcidChloride'}:peptide['sequence']=[value];index=0;site='N'
        peptide['modifications'].append({'site':site,'group':group})
    if kind=='paired_group_redistribution':
        peptide['sequence'].append(value);peptide['modifications'].append({'site':len(peptide['sequence']),'group':group})
    if kind=='site_bridge':
        for atom,partner in t['partner_atoms'].items():
            if t.get('external_partner'):break
            partner_value=(partner.get('residues') or ['G'])[0]
            if partner.get('position')=='start':continue
            peptide['sequence'].append(partner_value)
            if partner.get('group'):peptide['modifications'].append({'site':len(peptide['sequence']),'group':partner['group']})
    if r.get('bond_selector'):
        seq=['G','C','C','C','C','G'];bondkind=r['bond_selector']['kinds'][0]
        b={'kind':bondkind,'sites':[2,3]}
        if bondkind in {'thiosulfinate','thiosulfonate'}:b['oxidized_site']=2
        peptide={'sequence':seq,'nterm':'free','cterm':'acid','bonds':[b]};index=1
    present={'stage_tags':set(r.get('required_conditions',[]))|set(r.get('required_any_conditions',[])[:1]),
             'features':set(r.get('required_features',[])),'residues':set(peptide['sequence']),
             'acyl_donors':[],'carbonyl_donors':[]}
    fact=s.get('required_site_fact',r.get('required_site_fact'))
    if fact:
        for site in range(1,len(peptide['sequence'])+1):present['features'].add(f'{fact}:{site}')
        present['features'].add(f'{fact}:N')
    if kind=='capped_suffix':
        for i in range(1,len(peptide['sequence'])+1):present['features'].add(f"{t.get('boundary_fact','synthesis_boundary')}:{i}")
    if kind in {'o_acylate','n_acylate'}:present['acyl_donors']=[{'sequence':['A','G'],'nterm':'acetyl','cterm':'acid','modifications':[]}]
    if kind=='acyl_donor_self_condensation':
        peptide={'sequence':['G','P','Hyp'],'nterm':'free','cterm':'acid','modifications':[{'site':'N','group':'Fmoc'}]}
        present['acyl_donors']=[{'sequence':['G','P'],'nterm':'free','cterm':'acid','modifications':[{'site':'N','group':'Fmoc'}]}]
    if kind=='carbonyl_condensation':present['carbonyl_donors']=normalize_carbonyl_donors([{'name':'acetaldehyde','smiles':'CC=O'}])
    state=normalize_molecule(peptide)
    context=(present,{'stage_tags':set(),'features':set(),'residues':set()})
    return state,context

def independent_chain_formula(comp,peptide):
    f=merge({'H':2,'O':1},*(comp['residues'][x]['residue_formula'] for x in peptide['sequence']),comp['c_terminal_deltas'][peptide['cterm']],
            {'C':2,'H':2,'O':1} if peptide['nterm']=='acetyl' else {},
            *(comp['covalent_groups'][x['group']]['net_formula'] if 'group' in x else x['delta_formula'] for x in peptide['modifications']))
    return f

def oracle_delta(bundle,branch,state,event,context):
    t=branch['recipe']['transform'];kind=t['kind'];comp=bundle['components'];p=state['chains'][0]['peptide']
    if kind in {'hydrolytic_split'}:return {'H':2,'O':1}
    if kind=='prefix_dkp':return {}
    if kind=='scramble_disulfides':return {}
    if kind=='n_acetyl_change':return {} if t['move_to_oxygen'] else {'C':-2,'H':-2,'O':-1}
    if kind in {'o_acylate','n_acylate'}:return merge(independent_chain_formula(comp,context[0]['acyl_donors'][0]),{'H':-2,'O':-1})
    if kind=='acyl_donor_self_condensation':
        donor=independent_chain_formula(comp,context[0]['acyl_donors'][0])
        return merge(donor,donor,{'H':-2,'O':-1},negate(comp['covalent_groups']['Fmoc']['net_formula']),negate(independent_chain_formula(comp,p)))
    if kind=='carbonyl_condensation':return {'C':2,'H':2}
    if kind in {'site_bridge','trp_crosslink'}:
        delta=deepcopy(t['bond_delta'])
        if len(event['product_state']['chains'])>len(state['chains']):
            partner=t.get('external_partner') or p;delta=merge(delta,independent_chain_formula(comp,partner))
        return delta
    if 'delta_formula' in t:return t['delta_formula']
    selector=branch['recipe'].get('site_selector',{});site=next((int(x.rsplit(':',1)[-1]) for x in event['sites'] if x.rsplit(':',1)[-1].isdigit()),1)
    old=p['sequence'][site-1]
    group=selector.get('group') or (selector.get('groups') or [None])[0]
    if kind=='delete':
        count=t.get('count',1)
        removed=merge(*(comp['residues'][x]['residue_formula'] for x in p['sequence'][site-1:site-1+count]),
                      *(comp['covalent_groups'][x['group']]['net_formula'] if 'group' in x else x['delta_formula']
                        for x in p['modifications'] if type(x['site']) is int and site<=x['site']<site+count))
        return negate(removed)
    if kind=='duplicate_residue':return merge(*(comp['residues'][old]['residue_formula'] for _ in range(t.get('count',1))),
                                               comp['covalent_groups']['Phospho']['net_formula'] if t.get('insert_phosphorylated') else {})
    if kind=='substitute':
        new=t.get('residue') or t['residues_by_input'][old]
        return merge(comp['residues'][new]['residue_formula'],negate(comp['residues'][old]['residue_formula']),negate(comp['covalent_groups'][group]['net_formula']) if group else {})
    if kind=='replace_group':return merge(comp['covalent_groups'][t['product_group']]['net_formula'] if t.get('product_group') else t['replacement_formula'],negate(comp['covalent_groups'][group]['net_formula']))
    if kind=='capped_suffix':
        kept=event['product_state']['chains'][0]['origin_sites'];original=state['chains'][0]['origin_sites']
        cut=original.index(kept[0])
        return merge(t['cap_formula'],negate(merge(*(comp['residues'][x]['residue_formula'] for x in p['sequence'][:cut]))))
    raise ValueError('Missing independent oracle for '+kind)

def applicability_counterexample(state,branch):
    s=branch['recipe'].get('site_selector',{});bad=deepcopy(state)
    if branch['recipe'].get('bond_selector'):
        bad['bonds']=[];return bad
    if s.get('group') or s.get('groups'):
        for c in bad['chains']:
            c['peptide']['modifications']=[];c['peptide'].pop('monomer_state',None)
        return bad
    residues=set(s.get('residues') or s.get('motif') or [])
    chosen=next((x for x in ['G','A','P','L','S'] if x not in residues),None) if residues else (s.get('exclude_residues') or [None])[0]
    if chosen:
        for c in bad['chains']:
            c['peptide']['sequence']=[chosen]*len(c['peptide']['sequence']);c['peptide']['modifications']=[]
        return bad
    if s.get('cterm'):
        for c in bad['chains']:c['peptide']['cterm']='amide' if s['cterm']!='amide' else 'acid'
        return bad
    if s.get('nterm'):
        for c in bad['chains']:c['peptide']['nterm']='acetyl' if s['nterm']=='free' else 'free'
        return bad
    if s.get('next_backbone_amide'):
        for c in bad['chains']:c['peptide']['sequence']=['P']*len(c['peptide']['sequence'])
        return bad
    return None

def run():
    bundle,version=load_runtime();rows=[]
    for record in bundle['rules']:
        for branch in record['branches']:
            row={'rule_id':record['rule_id'],'branch_id':branch['branch_id'],'status':'fail'}
            try:
                state,ctx=fixture(bundle,branch)
                fn=bridge_transitions if branch['contract']['transform_family']=='bridge_edit' else lambda s,r,b,c:linear_transitions(bundle,s,r,b,c)
                events=list(fn(state,record,branch,ctx))
                if not events:raise ValueError('Positive fixture generated no event')
                for event in events:
                    expected=oracle_delta(bundle,branch,state,event,ctx)
                    if abs(m(expected)-m(event['delta_formula']))>1e-5:raise ValueError(f'Independent delta mismatch: {expected} != {event["delta_formula"]}')
                    if event['condition']=='conflict':raise ValueError('Positive fixture was excluded')
                req=branch['recipe'].get('required_conditions',[]);anyreq=branch['recipe'].get('required_any_conditions',[])
                reqfeatures=set(branch['recipe'].get('required_features',[]))
                fact=branch['recipe'].get('site_selector',{}).get('required_site_fact',branch['recipe'].get('required_site_fact'))
                if fact:reqfeatures|={x for x in ctx[0]['features'] if x.startswith(fact+':')}
                negative='not_required'
                if req or anyreq or reqfeatures:
                    bad=deepcopy(ctx);bad[0]['stage_tags']-=set(req+anyreq);bad[1]['stage_tags']=set(req+anyreq)
                    bad[0]['features']-=reqfeatures;bad[1]['features']=reqfeatures
                    neg=list(fn(state,record,branch,bad))
                    if not neg or any(e['condition']!='conflict' for e in neg):raise ValueError('Required-condition counterexample was not reported')
                    negative='pass'
                counter=applicability_counterexample(state,branch)
                app='unconstrained_substrate; positive applicability verified'
                if counter is not None:
                    if list(fn(counter,record,branch,ctx)):raise ValueError('Substrate/group/terminal counterexample still generated an event')
                    app='pass'
                row.update(status='pass',events=len(events),condition_negative=negative,applicability=app,oracle='Pyteomics + independent recipe arithmetic',fixture=state)
            except Exception as error:row['error']=str(error)
            rows.append(row)
    output={'library_hash':version,'scope':'practical compiled-branch mass/application and explicit-condition checks; exact chemistry/geometry is not claimed',
            'total':len(rows),'passed':sum(r['status']=='pass' for r in rows),'failed':sum(r['status']!='pass' for r in rows),'rows':rows}
    (ROOT/'rebuild/branch_acceptance.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='rows'},ensure_ascii=False))
    print(json.dumps([{'id':r['branch_id'],'error':r.get('error')} for r in rows if r['status']!='pass'],ensure_ascii=False))
    return output
if __name__=='__main__':run()
