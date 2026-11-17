"""Mass-family source leads, kept separate from executed structural candidates."""
import re
from copy import deepcopy
from observation_model import formula_mass
from structure_state import describe
from relative_case import applicability

def source_rows(bundle):
    rows=[]
    byid={r['rule_id']:r for r in bundle['rules']}
    def profile_for(ids,original=None):
        if original:return original
        recipes=[byid[i]['template'] for i in ids if i in byid and byid[i]['template']]
        if not recipes:return None
        recipe=recipes[0];selector=recipe.get('site_selector',{})
        profile={'stage_tags':recipe.get('required_conditions',[]),'stage_any':recipe.get('required_any_conditions',[]),
                 'features':recipe.get('required_features',[])}
        if selector.get('residues'):profile['any_residues']=selector['residues']
        return profile
    for record in bundle['rules']:
        r=record['legacy_record']
        if r is None:continue
        rows.append({'id':record['rule_id'],'name':r.get('name_cn',record['rule_id']),'formula':r.get('delta_formula'),
            'mass_status':r.get('mass_status'),'nominal':r.get('nominal_mass_shift_book'),
            'source':record.get('effective_source') or r.get('source') or r.get('source_basis'),
            'reference':r.get('mass_reference'),'evidence':r.get('certainty_notes') or r.get('diagnostic_evidence'),
            'profile':profile_for([record['rule_id']],r.get('relative_profile')),'kind':'rule_reference',
            'automatic':r.get('automatic_mass_event'),'origin':record['origin'],
            'branch_profiles':[{'stage_tags':b['recipe'].get('required_conditions',[]),
                'stage_any':b['recipe'].get('required_any_conditions',[]),'features':b['recipe'].get('required_features',[]),
                'any_residues':b['recipe'].get('site_selector',{}).get('residues',[]),
                'residues':b['recipe'].get('site_selector',{}).get('requires_residues',[])} for b in record['branches']]
                if len(record['branches'])>1 and not r.get('relative_profile') else []})
    for row in bundle['appendix_dispositions']:
        rows.append({'id':row['rule_id'],'name':row.get('connection') or row.get('review_note'),
            'formula':row.get('delta_formula'),'mass_status':'formula_reference','nominal':str(row.get('book_raw_value'))+' amu',
            'source':row['source'],'reference':row.get('reference_state'),'evidence':row.get('evidence_level'),
            'kind':row.get('record_kind','appendix_reference'),'origin':'book_appendix','profile':profile_for(row.get('legacy_rule_ids',[])),
            'intrinsic_charge':row.get('intrinsic_charge',0)})
    return rows

def source_applicability(row,context):
    profiles=row.get('branch_profiles') or [row.get('profile')]
    fits=[applicability(row['id'],context,p) for p in profiles]
    ranks={'context_supported':3,'context_partial':2,'unreviewed_context':1,'context_conflict':0}
    return max(fits,key=lambda x:ranks[x[0]])

def lookup_delta(bundle,value,tolerance=None,sequence=None,tags=(),absent=()):
    ctx=({'stage_tags':set(tags),'features':set(),'residues':set(sequence or [])},
         {'stage_tags':set(absent),'features':set(),'residues':set(bundle['components']['residues'])-set(sequence) if sequence else set()})
    found=[]
    for row in source_rows(bundle):
        nominal=re.match(r'^\s*([+−–-]?\d+)\s*amu',str(row['nominal'] or ''))
        theoretical=formula_mass(row['formula']) if isinstance(row['formula'],dict) and not row.get('intrinsic_charge') else None
        if tolerance is None:
            matched=(int(nominal[1].replace('−','-').replace('–','-'))==round(value)) if nominal else theoretical is not None and round(theoretical)==round(value)
            basis='book_integer_source_lookup_without_measurement_tolerance' if nominal else 'rounded_formula_source_hint_without_measurement_tolerance'
        else:
            matched=theoretical is not None and abs(value-theoretical)<=tolerance
            basis='formula_delta_source_lookup'
        if not matched:continue
        fit,yes,missing,no=source_applicability(row,ctx)
        found.append({**row,'theoretical_delta':theoretical,'match_basis':basis,'context_status':fit,
                      'missing_conditions':missing,'conflicting_conditions':no,'classification':'source_lead_only',
                      'reference_status':'unverified_reference','structure_confirmed':False})
    found.sort(key=lambda x:(x['context_status']=='context_conflict',x['context_status'] in {'context_partial','unreviewed_context'},x['id']))
    return {'lookup_mode':True,'value':value,'tolerance':tolerance,'source_leads':found,'state_search_performed':False,
            'precision_note':'仅按整数偏差粗查（原书值或净式舍入），未判断高分辨吻合' if tolerance is None else '按指定误差查净质量家族，参照和结构仍需核对'}

def mass_source_leads(bundle,observation,options,references,context,target=None):
    ctx=deepcopy(context)
    if target:
        ctx[0]['residues']={r for c in target['chains'] for r in c['peptide']['sequence']}
        ctx[1]['residues']=set(bundle['components']['residues'])-ctx[0]['residues']
    results={}
    for row in source_rows(bundle):
        if row['mass_status']=='not_product_mass_rule':continue
        for option in options:
            f=row['formula'];basis=option['basis'];delta=formula_mass(f,basis) if isinstance(f,dict) else 0 if row['mass_status']=='isomeric' else None
            if row.get('intrinsic_charge') or row['id']=='SUP-039':delta=None
            for reference in references:
                refmass=formula_mass(reference['formula'],basis)
                if option['kind']=='absolute':
                    coordinate=option['coordinate']
                    if coordinate.intrinsic_charge:continue
                    measured=coordinate.elemental_reference()-refmass;tol=coordinate.elemental_tolerance()
                else:measured=option['value'];tol=option['tolerance']
                evidence='formula_mass_reference'
                if delta is not None:match=abs(measured-delta)<=tol;error=measured-delta
                else:
                    nominal=re.match(r'^\s*([+−–-]?\d+)\s*amu',str(row['nominal'] or ''))
                    match=observation.get('kind')=='delta_approx' and nominal is not None and abs(measured-int(nominal[1].replace('−','-').replace('–','-')))<=tol
                    error=None;evidence='book_integer_clue_only'
                if not match:continue
                fit,yes,missing,no=source_applicability(row,ctx)
                key=(row['id'],option['label'],tuple(reference['chain_ids']))
                results[key]={**row,'theoretical_delta':delta,'error_da':error,'mass_basis':basis,
                    'match_basis':evidence,'reference_status':'source_reference_not_asserted_equal_to_current_target',
                    'context_status':fit,'supported_conditions':yes,'missing_conditions':missing,'conflicting_conditions':no,
                    'ion_assumption':option['label'],'reference_chain_ids':reference['chain_ids'],
                    'classification':'source_lead_only','structure_confirmed':False}
    return sorted(results.values(),key=lambda x:x['id'])
