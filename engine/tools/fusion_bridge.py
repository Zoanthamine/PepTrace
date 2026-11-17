"""Pinned DSH arithmetic/source review with the canonical state core as executor."""
from copy import deepcopy
from functools import lru_cache
import hashlib,json,sys
from pathlib import Path
from knowledge_source import ROOT,load_runtime
from state_diagnosis import diagnose_state_case,observation_options,reference_components,prepare_context
from structure_state import normalize_molecule
from observation_model import formula_mass
from state_transformations import conditions
sys.path.insert(0,str(ROOT))
from fusion_vendor.dsh_core.components import load_library
from fusion_vendor.dsh_core.rules import load_rules
from fusion_vendor.dsh_core.formula import Formula
from fusion_vendor.dsh_core.matching import match_delta

@lru_cache(maxsize=1)
def peer_assets():
    provenance=json.loads((ROOT/'fusion_vendor/provenance.json').read_text(encoding='utf-8'))
    for name,expected in provenance['files_sha256'].items():
        if hashlib.sha256((ROOT/'fusion_vendor/dsh_core'/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Pinned DSH module changed; revalidate release: '+name)
    lib=load_library(str(ROOT/'多肽质量组件.json'));rules=load_rules(lib,kb_root=ROOT)
    if rules.problems:raise ValueError('DSH independent component/rule arithmetic failed: '+str(rules.problems))
    return lib,rules

def diagnose_fused_case(case):
    result=diagnose_state_case(case);bundle,version=load_runtime();lib,peer=peer_assets()
    target=normalize_molecule(case['intended_product']);starting=normalize_molecule(case['starting_material']) if case.get('starting_material') else None
    ctx=prepare_context(case,target,starting);byid={r['rule_id']:r for r in bundle['rules']}
    all_reviews=[];mass_checks=[]
    for observation in result['observations']:
        raw=observation['input'];options=observation_options(raw);refs=reference_components(case,raw,target,starting)
        found={}
        for option in options:
            if option['basis']!='monoisotopic':continue
            for ref in refs:
                if option['kind']=='absolute':
                    coordinate=option['coordinate']
                    if coordinate.intrinsic_charge or ref['intrinsic_charge']:continue
                    delta=coordinate.elemental_reference()-formula_mass(ref['formula']);tol=coordinate.elemental_tolerance()
                else:
                    if option['kind']=='pair' and (option['first'].intrinsic_charge or option['second'].intrinsic_charge):continue
                    delta=option['value'];tol=option['tolerance']
                for match in match_delta(peer.mass_matchable(),delta,tol):
                    rule=match.rule
                    if rule.mass_status=='not_product_mass_rule' or rule.polymer_relation or rule.rule_id=='SUP-039':continue
                    if match.unit_ratio>case.get('max_events',2):continue
                    if match.grade.startswith('nominal'):
                        if raw['kind']!='delta_approx':continue
                    elif abs(match.error_da)>tol:continue
                    row=byid[rule.rule_id];fits=[{'branch_id':b['branch_id'],**conditions(b['recipe'],ctx)} for b in row['branches']]
                    found[(rule.rule_id,option['label'],tuple(ref['chain_ids']))]={
                        'id':rule.rule_id,'name':rule.name_cn,'theoretical_delta':match.predicted,'error_da':match.error_da,
                        'mass_grade':match.grade,'unit_ratio':match.unit_ratio,'source':row['effective_source'],
                        'reference':rule.mass_reference,'diagnostic_evidence':rule.diagnostic_evidence,
                        'competing_explanations':rule.competing,'mitigation':rule.mitigation,
                        'canonical_branch_conditions':fits,'classification':'peer_source_lead_only',
                        'structure_confirmed':False,'origin':rule.origin,'ion_assumption':option['label']}
        reviews=list(found.values());observation['dsh_source_review']=reviews
        core_ids={x['id'] for x in observation['source_leads']};peer_ids={x['id'] for x in reviews}
        all_reviews.append({'core_source_ids':sorted(core_ids),'peer_source_ids':sorted(peer_ids),
            'peer_only_source_ids':sorted(peer_ids-core_ids),'core_only_source_ids':sorted(core_ids-peer_ids),
            'interpretation':'Source lookup differences retained for review; additive peer matches never become executed structural candidates.'})
        for h in observation['candidates']+observation['condition_conflicts']:
            if h['mass_basis']=='monoisotopic':
                independent=Formula.of(h['component']['formula']).mass()
                error=independent-h['component']['element_mass'];mass_checks.append(abs(error)<1e-5)
                if abs(error)>=1e-5:raise ValueError('DSH independent whole-product formula check failed')
            notes=[]
            for rid in sorted({e['rule_id'] for p in h['paths'] for e in p}):
                row=byid[rid];legacy=row['legacy_record'] or {}
                notes.append({'rule_id':rid,'name':legacy.get('name_cn',rid),'origin':row['origin'],
                    'source':row['effective_source'],'diagnostic_evidence':legacy.get('diagnostic_evidence',[]),
                    'competing_explanations':legacy.get('competing_explanations',[]),'mitigation':legacy.get('mitigation',[]),
                    'branch_evidence':[b['recipe'].get('diagnostic_evidence',[]) for b in row['branches']]})
            h['source_notes']=notes
    result['engine']='codex-dsh-fused-1'
    result['fusion_audit']={'canonical_library_hash':version,'dsh_component_path':'多肽质量组件.json',
        'dsh_rule_root':'.','dsh_rules':len(peer),'canonical_rules':len(bundle['rules']),
        'dsh_data_problems':peer.problems,'whole_product_formula_checks':len(mass_checks),
        'whole_product_formula_checks_pass':all(mass_checks),'source_comparisons':all_reviews,
        'formal_executor':'Codex state core','condition_authority':'canonical branch ALL/ANY predicates',
        'not_used_as_executor':['DSH legacy additive combination search','DSH old conditions_map','DSH implicit observation defaults']}
    return result

def adapt_dsh_case(raw):
    """Support explicit legacy JSON; reject defaults that would hide missing input."""
    if set(raw)-{'case_id','sequence','nterm','cterm','modifications','conditions','observed','notes','fixed_charge'}:
        raise ValueError('Unknown legacy case fields')
    for key in ('sequence','nterm','cterm','observed'):
        if key not in raw:raise ValueError('Legacy adapter requires explicit '+key)
    if raw.get('fixed_charge',0):raise ValueError('Declare fixed-charge atoms with the native state schema')
    obs=raw['observed'];kind=obs.get('kind')
    if set(obs)-{'kind','value','z','z2','adduct','tol_da','tol_ppm','value2','note'}:raise ValueError('Unknown legacy observation fields')
    if not any(x in obs for x in ('tol_da','tol_ppm')):raise ValueError('Legacy adapter requires an explicit tolerance')
    if kind=='delta':
        if 'tol_da' not in obs:raise ValueError('Relative observation requires a Da window')
        converted={'kind':'delta_mono','value':obs['value'],'tolerance':obs['tol_da'],'tolerance_unit':'Da'}
    elif kind in {'mz','neutral_mass'}:
        converted={'kind':'mz' if kind=='mz' else 'neutral_mono','value':obs['value'],
            'tolerance':obs.get('tol_da',obs.get('tol_ppm')),'tolerance_unit':'Da' if 'tol_da' in obs else 'ppm'}
        if kind=='mz':
            if 'z' not in obs or 'adduct' not in obs:raise ValueError('Legacy mz requires explicit z and adduct')
            z=obs['z'];adduct=obs['adduct']
            if adduct not in {'H','Na','K','-H'}:raise ValueError('Use native ion schema for this adduct')
            converted.update(charge=z,ions={adduct:abs(z) if adduct in {'H','-H'} else 1})
    else:raise ValueError('Use native peak_pair schema for explicit dual-peak uncertainties')
    from natural_case import ALIASES
    cond=raw.get('conditions',{})
    if set(cond)-{'reported_reagents','denied_reagents','tags','stage','notes'}:raise ValueError('Unknown legacy condition fields')
    if cond.get('tags') or cond.get('stage'):raise ValueError('Map legacy tags/stage explicitly with the native reaction schema')
    def tags(values):return sorted({ALIASES.get(str(v).lower(),str(v).lower()) for v in values})
    reaction={'description':raw.get('notes') or cond.get('notes') or 'Legacy DSH structured input',
        'stage_tags':tags(cond.get('reported_reagents',[])),'stage_tags_absent':tags(cond.get('denied_reagents',[]))}
    return {'case_id':raw.get('case_id','legacy-dsh'),'intended_product':{k:deepcopy(raw[k]) for k in ('sequence','nterm','cterm','modifications') if k in raw},
        'reaction':reaction,'observations':[converted]}

