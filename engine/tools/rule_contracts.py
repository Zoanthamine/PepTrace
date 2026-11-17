"""Compile explicit Codex recipes into state/condition/consumption contracts."""
from copy import deepcopy

FAMILIES={
 'sequence_edit':{'delete','duplicate_residue','insert_before','insert_gly_before','substitute','capped_suffix','prefix_dkp','acyl_donor_self_condensation'},
 'group_edit':{'local_delta','n_terminal_delta','c_terminal_delta','remove_group','replace_group','oxidize_protected_thioether','n_acetyl_change'},
 'connection_edit':{'site_bridge','trp_crosslink','o_acylate','n_acylate','backbone_amide_to_ester','hydrolytic_split','migrate_group_to_previous','transfer_group_to_amine','transfer_group_to_oxygen','carbonyl_condensation'},
 'terminal_edit':{'c_terminal_amide_substitution','c_terminal_pyrrolidinamide','c_terminal_acid'},
 'bridge_edit':{'convert_sulfur_bridge','oxidize_sulfur_bridge','reduce_sulfur_bridge','hydrolyze_disulfide','hydrolyze_thiosulfinate','oxidative_bridge_cleavage','disulfide_beta_elimination','disulfide_alpha_elimination','scramble_disulfides'},
 'state_edit':{'epimerize_alpha','paired_group_redistribution'}}

def compile_contract(recipe):
    transform=recipe['transform'];kind=transform['kind']
    family=next((name for name,kinds in FAMILIES.items() if kind in kinds),None)
    if family is None:raise ValueError('Unclassified transformation '+kind)
    selector=recipe.get('site_selector',recipe.get('bond_selector',{}))
    groups=[selector['group']] if 'group' in selector else selector.get('groups',[])
    consumed=[]
    if groups:consumed.append({'resource':'bound_covalent_group','groups':groups,'selector':deepcopy(selector)})
    if family=='bridge_edit':consumed.append({'resource':'selected_bridge_state','kinds':selector.get('kinds',[])})
    elif kind in {'site_bridge','trp_crosslink','o_acylate','n_acylate'}:
        consumed.append({'resource':'reacting_connection_endpoints','source_atom':transform.get('source_atom'),
                         'partner_atoms':deepcopy(transform.get('partner_atoms',{}))})
    elif kind=='epimerize_alpha':consumed.append({'resource':'selected_alpha_configuration','update':'relative_inversion_or_retention'})
    else:consumed.append({'resource':'selected_input_structural_state','selector':deepcopy(selector)})
    prerequisites=deepcopy(recipe.get('required_features',[]))
    if recipe.get('required_site_fact') or selector.get('required_site_fact'):
        prerequisites.append({'site_fact':selector.get('required_site_fact',recipe.get('required_site_fact')),
                              'site':selector.get('fact_site','selected_site')})
    return {'schema_version':'state-rule-contract-1','transform_family':family,
            'reference_state':recipe['reference_state'],'substrate_predicates':deepcopy(selector),
            'conditions':{'required_all':recipe.get('required_conditions',[]),'required_any':recipe.get('required_any_conditions',[]),
                          'promoting_any':recipe.get('supporting_any_conditions',[]),'unknown_policy':'retain','conflict_policy':'report'},
            'prerequisites':prerequisites,'consumed_resources':consumed,
            'exclusive_policy':'current_state_must_supply_each_consumed_resource',
            'repeat_limit':{'same_input_binding':1,'other_sites':'enumerate','new_product_state':'recheck_predicates'},
            'product_state_operator':deepcopy(transform),
            'intrinsic_charge_change':transform.get('intrinsic_charge',0),
            'source_refs':deepcopy(recipe['source_refs']),'evidence_level':recipe['evidence_level']}

def validate_contract(contract):
    if contract['conditions']['unknown_policy']!='retain':raise ValueError('Unknown conditions must be retained')
    if not contract['consumed_resources']:raise ValueError('A reaction must declare structural resource consumption')
    if contract['repeat_limit']['same_input_binding']!=1:raise ValueError('Duplicate replay at identical input binding')
    if contract['transform_family'] not in FAMILIES:raise ValueError('Unknown transformation family')
    if type(contract['intrinsic_charge_change']) is not int:raise ValueError('Intrinsic charge must be an integer')

def legacy_route(record):
    if record['template']:return {'kind':'executable_template','branch_ids':[b['branch_id'] for b in record['branches']]}
    disposition=record.get('disposition') or {}
    if disposition.get('status')=='compiled_fact':return {'kind':'source_fact','reason':disposition.get('unresolved_reason')}
    return {'kind':'source_reference_only','reason':disposition.get('unresolved_reason') or
            'Source record has no compiled structural recipe; do not turn a lookup delta into a reaction.'}
