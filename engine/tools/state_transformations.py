"""Apply compiled recipes to current molecular states, preserving connections and residue identities."""
from copy import deepcopy
from case_input import evaluate_requirements
from linear_rules import instantiate, selected_sites
from atom_model import query_atoms,consumed_atoms,atom_options,mod_atoms,site_matches,default_mod_atom
from structure_state import (add, chain_formula, components_of, encoded, normalize_molecule,
                             projection_formula, state_formula, validate_state)

class UnsupportedTransition(ValueError):pass

def local_context(context,chain):
    out=deepcopy(context);cid=chain['id'];origins=chain['origin_sites']
    for family in ('features',):
        for part in out:
            result=set()
            for feature in part[family]:
                if '|' in feature:
                    prefix,feature=feature.split('|',1)
                    if prefix!=cid:continue
                head,sep,tail=feature.rpartition(':')
                if sep and tail.isdigit():
                    original=f'{cid}:{tail}'
                    if original in origins:result.add(head+':'+str(origins.index(original)+1))
                else:result.add(feature)
            part[family]=result
    return out

def site_origin(chain,site):
    return chain['origin_sites'][site-1] if type(site) is int else chain['id']+':'+str(site)

def conditions(recipe,context):
    support,missing,conflict=evaluate_requirements(recipe.get('required_conditions',[]),recipe.get('required_any_conditions',[]),
                                                 context[0]['stage_tags'],context[1]['stage_tags'])
    yes,unknown,no=evaluate_requirements(recipe.get('required_features',[]),[],context[0]['features'],context[1]['features'])
    support+=yes;missing+=unknown;conflict+=no
    support+=sorted(set(recipe.get('supporting_any_conditions',[]))&context[0]['stage_tags'])
    return {'condition':'conflict' if conflict else 'unknown' if missing else 'compatible',
            'supported_conditions':support,'missing_conditions':missing,'conflicting_conditions':conflict}

def pack_for(bundle,record,branch):
    recipe={k:v for k,v in branch['recipe'].items() if k!='alternate_branches'}
    recipe['_branch_id']=branch['branch_id']
    return {'schema_version':1,'source':bundle['pack_metadata'][record['compatibility_pack']]['source'],'rules':[recipe]}

def transfer_bonds(state,old_chain,new_chains):
    lookup={origin:(chain['id'],i) for chain in new_chains for i,origin in enumerate(chain['origin_sites'],1)}
    bonds=deepcopy(state['bonds'])
    for bond in bonds:
        endpoints=bond['endpoints']+([bond['oxidized_endpoint']] if 'oxidized_endpoint' in bond else [])
        for endpoint in endpoints:
            if endpoint['chain']!=old_chain['id']:continue
            site=endpoint['site']
            if type(site) is int:
                origin=old_chain['origin_sites'][site-1]
                if origin not in lookup:raise UnsupportedTransition('Sequence edit would remove a connected atom without a defined bond-cleavage product')
                endpoint['chain'],endpoint['site']=lookup[origin]
            elif len(new_chains)==1:endpoint['chain']=new_chains[0]['id']
            elif site=='N':endpoint['chain']=next(c['id'] for c in new_chains if old_chain['origin_sites'][0] in c['origin_sites'])
            else:endpoint['chain']=next(c['id'] for c in new_chains if old_chain['origin_sites'][-1] in c['origin_sites'])
    return bonds

def chain_product(old_chain,product,cid,recipe,event):
    allowed={'sequence','nterm','cterm','modifications','alpha_stereochemistry','monomer_state'}
    peptide={k:deepcopy(v) for k,v in product.items() if k in allowed}
    for modification in peptide['modifications']:
        if 'delta_formula' in modification and 'atom' not in modification and 'atoms' not in modification:
            modification['atom']=default_mod_atom(peptide,modification)
    origins=product.get('_origin_sites')
    if origins is None:
        if len(peptide['sequence'])==len(old_chain['origin_sites']):origins=old_chain['origin_sites']
        else:origins=[f'{cid}:created:{i}' for i in range(1,len(peptide['sequence'])+1)]
    if recipe['transform']['kind']!='epimerize_alpha':
        old_stereo=old_chain['peptide'].get('alpha_stereochemistry',{})
        preserved={str(i+1):old_stereo[str(old_chain['origin_sites'].index(origin)+1)]
                   for i,origin in enumerate(origins) if origin in old_chain['origin_sites'] and
                   str(old_chain['origin_sites'].index(origin)+1) in old_stereo and
                   peptide['sequence'][i] not in {'G','Aib','beta-Ala','Aoa','Dha','MeDha'}}
        if preserved:peptide['alpha_stereochemistry']=preserved
        else:peptide.pop('alpha_stereochemistry',None)
    unknown_charge=old_chain.get('unlocalized_intrinsic_charge',old_chain['intrinsic_charge'] if not old_chain.get('charge_locations') else 0)
    if unknown_charge and origins!=old_chain['origin_sites']:
        raise UnsupportedTransition('Charged fragmentation requires the locations of permanent charged groups')
    charge_locations=[deepcopy(x) for x in old_chain.get('charge_locations',[]) if x['origin'] in origins]
    for location in charge_locations:
        if location['atom']=='backbone_N' and origins.index(location['origin'])==0:location['atom']='Nalpha'
    delta_charge=recipe['transform'].get('intrinsic_charge',0)
    if delta_charge:
        selected=event['site'];point=selected if type(selected) is int or selected in ('N','C') else next((s for s in event['consumed'] if type(s) is int),None)
        atom_set=query_atoms(recipe,old_chain['peptide'],point) if point is not None else set()
        if point is not None and len(atom_set)==1:
            origin=old_chain['origin_sites'][0 if point=='N' else -1 if point=='C' else point-1]
            charge_locations.append({'origin':origin,'atom':next(iter(atom_set)),'charge':delta_charge,'source_rule':recipe['rule_id']})
        else:unknown_charge+=delta_charge
    annotations=deepcopy(old_chain['annotations'])
    for annotation in annotations:
        if 'consumed_sites' in annotation:
            origin_refs=annotation.get('consumed_origins',[site_origin(old_chain,s) for s in annotation['consumed_sites']])
            annotation['consumed_origins']=origin_refs
            annotation['consumed_sites']=[origins.index(o)+1 if o in origins else o.rsplit(':',1)[-1]
                                          for o in origin_refs if o in origins or o.endswith(':N') or o.endswith(':C')]
        if annotation.get('kind')=='product_provenance':
            original=annotation.get('site_origin')
            if original in origins:annotation['site']=origins.index(original)+1
    if 'cycle' in product or 'bond_edit' in product or 'linkage_edit' in product:
        annotations.append({'branch_id':recipe['_branch_id'],'rule_id':recipe['rule_id'],
                            'kind':recipe['transform']['kind'],'cycle':product.get('cycle'),
                            'bond_edit':product.get('bond_edit'),'linkage_edit':product.get('linkage_edit'),
                            'consumed_sites':[x for x in event['consumed'] if type(x) is int or x in ('N','C')],
                            'consumed_origins':[site_origin(old_chain,x) for x in event['consumed'] if type(x) is int or x in ('N','C')]})
        annotations[-1]['consumed_atoms']=consumed_atoms(recipe,old_chain['peptide'],event)
        if recipe['transform'].get('consume_nterm'):
            annotations[-1]['consumed_sites'].append('N');annotations[-1]['consumed_origins'].append(old_chain['id']+':N')
    if recipe['transform'].get('unresolved_linkage'):
        annotations.append({'kind':'unresolved_linkage','site':event['site'],'description':recipe['transform']['unresolved_linkage']})
    # Store provenance for defined subsequent state conversions, separately from the formula.
    if recipe.get('_keep_provenance',True):annotations.append({'branch_id':recipe['_branch_id'],'rule_id':recipe['rule_id'],
                        'kind':'product_provenance','site':event['site'],
                        'modification':recipe['transform'].get('modification'),
                        'site_origin':site_origin(old_chain,event['site']) if type(event['site']) is int or event['site'] in ('N','C') else None})
    relative=deepcopy(old_chain.get('relative_alpha_stereo',{}))
    if recipe['transform']['kind']=='epimerize_alpha' and recipe['transform'].get('invert',True):
        site=event['site']
        if type(site) is int and str(site) not in old_chain['peptide'].get('alpha_stereochemistry',{}):
            origin=site_origin(old_chain,site)
            if origin in relative:relative.pop(origin)
            else:relative[origin]='inverted_relative_to_input'
    return {'id':cid,'peptide':peptide,'origin_sites':deepcopy(origins),
            'intrinsic_charge':sum(x['charge'] for x in charge_locations)+unknown_charge,
            'annotations':annotations,**({'relative_alpha_stereo':relative} if relative else {}),
            **({'charge_locations':charge_locations} if charge_locations else {}),
            **({'unlocalized_intrinsic_charge':unknown_charge} if unknown_charge else {})}

def replace_chain(state,chain,products,recipe,events):
    output=deepcopy(state);new=[]
    for i,(product,event) in enumerate(zip(products,events),1):
        cid=chain['id'] if len(products)==1 else chain['id']+f'.fragment{i}'
        new.append(chain_product(chain,product,cid,recipe,event))
    output['chains']=[c for c in output['chains'] if c['id']!=chain['id']]+new
    output['bonds']=transfer_bonds(state,chain,new)
    for new_chain,product in zip(new,products):
        if product.get('topology')=='head_tail_cycle':
            output['bonds'].append({'kind':'amide','endpoints':[{'chain':new_chain['id'],'site':'N','atom':'Nalpha'},
                                                             {'chain':new_chain['id'],'site':'C','atom':'carbonyl_C'}],
                                    'delta_formula':{'H':-2,'O':-1}})
    output['chains'].sort(key=lambda c:c['id']);output['bonds'].sort(key=encoded)
    validate_state(output)
    return output

def continuation_projections(bundle,record,chain,target):
    paths=[p for p in bundle['continuation_paths']['paths'] if p['next_rule_id']==record['rule_id']]
    for path in paths:
        for mod in target['modifications']:
            site=mod['site']
            provenance=[a for a in chain['annotations'] if a.get('kind')=='product_provenance'
                        and (a.get('site')==site or isinstance(a.get('site'),list) and site in a['site'])]
            if mod.get('group') not in path.get('previous_groups',[]) and not any(a['rule_id'] in path.get('previous_rule_ids',[]) and
                    a.get('modification')==mod.get('name') for a in provenance):continue
            projection=deepcopy(target)
            projection['modifications']=[m for m in projection['modifications'] if m['site']!=site]
            yield projection,site,path

def linear_transitions(bundle,state,record,branch,context):
    recipe={k:v for k,v in branch['recipe'].items() if k!='alternate_branches'};recipe['_branch_id']=branch['branch_id']
    recipe['_keep_provenance']=any(record['rule_id'] in p.get('previous_rule_ids',[]) for p in bundle['continuation_paths']['paths'])
    kind=recipe['transform']['kind']
    if kind in {'site_bridge','trp_crosslink','o_acylate','n_acylate'}:
        yield from connection_transitions(state,record,branch,context)
        return
    for chain in state['chains']:
        p=chain['peptide'];s=recipe['site_selector'];codes=set(p['sequence'])
        if s.get('residues') and not codes&set(s['residues']):continue
        if s.get('requires_residues') and not set(s['requires_residues'])<=codes:continue
        if s.get('exclude_residues') and codes<=set(s['exclude_residues']):continue
        if s.get('sequence_length') and len(p['sequence'])!=s['sequence_length']:continue
        if s.get('nterm') and p['nterm']!=s['nterm'] or s.get('cterm') and p['cterm']!=s['cterm']:continue
        groups={m.get('group') for m in p['modifications']}
        if s.get('group') and s['group'] not in groups:continue
        if s.get('groups') and not set(s['groups'])&groups:continue
        target=deepcopy(chain['peptide']);target['_origin_sites']=chain['origin_sites']
        target['_atom_model']=True
        ctx=local_context(context,chain)
        pack=pack_for(bundle,record,branch)
        projections=[(target,None,None)]+list(continuation_projections(bundle,record,chain,target))
        for projected,continuation_site,path in projections:
            def formula_of(p):
                if p is not projected:return p['sequence'],chain_formula(p),set()
                if continuation_site is None:return projection_formula(state,chain['id'],p,recipe)
                clean_chain=deepcopy(chain)
                clean_chain['annotations']=[a for a in clean_chain['annotations'] if not (a.get('rule_id') in path.get('previous_rule_ids',[]))]
                clean_state={**state,'chains':[clean_chain if c['id']==chain['id'] else c for c in state['chains']]}
                return projection_formula(clean_state,chain['id'],p,recipe)
            generated=list(instantiate(pack,projected,ctx,formula_of,add))
            if continuation_site is not None:
                generated=[e for e in generated if e['site']==continuation_site]
            if kind in {'hydrolytic_split','prefix_dkp'}:
                groups={}
                for event in generated:groups.setdefault(encoded(event['consumed']),[]).append(event)
                for paired in groups.values():
                    if len(paired)!=2:raise UnsupportedTransition('Fragment recipe lacks both coproduced fragments')
                    output=replace_chain(state,chain,[e['product_structure'] for e in paired],recipe,paired)
                    yield make_transition(state,output,chain,record,branch,paired[0],paired[0]['consumed'],
                                          extra={'coproducts':[e['product'] for e in paired]})
                continue
            for event in generated:
                output=replace_chain(state,chain,[event['product_structure']],recipe,[event])
                if continuation_site is not None:
                    oc=next(c for c in output['chains'] if c['id']==chain['id'])
                    # The product conversion replaces the previous structural ring/state, rather than accumulating two rings.
                    oc['annotations']=[a for a in oc['annotations'] if a.get('rule_id') not in path.get('previous_rule_ids',[])]
                yield make_transition(state,output,chain,record,branch,event,event['consumed'],
                                      extra={'continuation_source':path} if path else {})

def make_transition(before,after,chain,record,branch,event,sites,extra=None):
    delta=add(state_formula(after),{e:-n for e,n in state_formula(before).items()})
    resources=[site_origin(chain,s) if type(s) is int or s in ('N','C') else encoded(s) for s in sites]
    return {'rule_id':record['rule_id'],'branch_id':branch['branch_id'],'sites':resources,'product_state':after,
            'delta_formula':delta,'contract':branch['contract'],
            'condition':event['condition'],'supported_conditions':event.get('supported_conditions',[]),
            'missing_conditions':event.get('missing_conditions',[]),'conflicting_conditions':event.get('conflicting_conditions',[]),
            'source':event.get('source',record.get('effective_source')),'reference_state':event.get('reference_state'),
            'product':event.get('product',branch['recipe']['product_description']),
            'evidence_level':branch['recipe']['evidence_level'],**(extra or {})}

def bridge_transitions(state,record,branch,context):
    recipe=branch['recipe'];selector=recipe['bond_selector'];transform=recipe['transform'];kind=transform['kind']
    if kind=='scramble_disulfides':
        yield from network_transitions(state,record,branch,context)
        return
    for i,bond in enumerate(state['bonds']):
        if bond['kind'] not in selector['kinds']:continue
        ends=bond['endpoints']
        choices=ends if kind in {'disulfide_beta_elimination','disulfide_alpha_elimination','hydrolyze_disulfide'} or kind=='oxidize_sulfur_bridge' and bond['kind']=='disulfide' else [bond.get('oxidized_endpoint',ends[0])]
        if kind=='convert_sulfur_bridge' and transform['product_kind']=='lanthionine':choices=ends
        for reacting in choices:
            output=deepcopy(state);b=output['bonds'][i]
            def group(endpoint,name):
                c=next(c for c in output['chains'] if c['id']==endpoint['chain'])
                c['peptide']['modifications'].append({'site':endpoint['site'],'group':name})
            if kind=='oxidize_sulfur_bridge':
                b.update(kind=transform['product_kind'],oxidized_endpoint=deepcopy(reacting))
                b['delta_formula']={'H':-2,'O':1 if b['kind']=='thiosulfinate' else 2}
            elif kind=='convert_sulfur_bridge':
                from structure_state import CONNECTION_DELTAS
                b.update(kind=transform['product_kind'],delta_formula=deepcopy(CONNECTION_DELTAS[transform['product_kind']]))
                b.pop('oxidized_endpoint',None)
                if b['kind']=='lanthionine':b['lost_sulfur_endpoint']=deepcopy(reacting)
            elif kind=='reduce_sulfur_bridge':output['bonds'].pop(i)
            elif kind=='oxidative_bridge_cleavage':
                output['bonds'].pop(i)
                for e in ends:group(e,'CysSulfonic')
            elif kind=='hydrolyze_thiosulfinate':
                output['bonds'].pop(i)
                for e in ends:group(e,'CysSulfenic')
            elif kind in {'disulfide_alpha_elimination','hydrolyze_disulfide'}:
                output['bonds'].pop(i);group(reacting,'CysThioaldehyde' if kind=='disulfide_alpha_elimination' else 'CysSulfenic')
            elif kind=='disulfide_beta_elimination':
                output['bonds'].pop(i)
                c=next(c for c in output['chains'] if c['id']==reacting['chain'])
                c['peptide']['sequence'][reacting['site']-1]='Dha'
                c['peptide'].get('alpha_stereochemistry',{}).pop(str(reacting['site']),None)
                group(next(e for e in ends if e!=reacting),'Persulfide')
            else:raise UnsupportedTransition('Unimplemented global bridge transform '+kind)
            validate_state(output)
            delta=add(state_formula(output),{e:-n for e,n in state_formula(state).items()})
            if 'delta_formula' in transform and delta!=transform['delta_formula']:
                raise ValueError('Global bridge element balance disagrees with explicit recipe '+branch['branch_id'])
            fit=conditions(recipe,context)
            yield {'rule_id':record['rule_id'],'branch_id':branch['branch_id'],
                   'sites':[e['chain']+':'+str(e['site']) for e in ends],
                   'reacting_endpoint':reacting,'product_state':output,'delta_formula':delta,'contract':branch['contract'],
                   'source':recipe['source_refs'],'reference_state':recipe['reference_state'],
                   'product':recipe['product_description'],'evidence_level':recipe['evidence_level'],**fit}

def transitions(bundle,state,context,rule_ids=None):
    unsupported=[];events=[]
    for record in bundle['rules']:
        if rule_ids is not None and record['rule_id'] not in rule_ids:continue
        for branch in record['branches']:
            try:
                if branch['contract']['transform_family']=='bridge_edit':events.extend(bridge_transitions(state,record,branch,context))
                else:events.extend(linear_transitions(bundle,state,record,branch,context))
            except UnsupportedTransition as error:
                unsupported.append({'rule_id':record['rule_id'],'branch_id':branch['branch_id'],'reason':str(error)})
    return events,unsupported

def endpoint_available(state,chain,site,atom,allowed_groups=()):
    p=chain['peptide']
    if site=='N' and p['nterm']!='free' or site=='C' and p['cterm']!='acid':return False
    residue=p['sequence'][0 if site=='N' else -1 if site=='C' else site-1];wanted=atom_options(atom,residue,site);length=len(p['sequence'])
    if any(site_matches(m['site'],site,length) and wanted&mod_atoms(p,m) and m.get('group') not in allowed_groups for m in p['modifications']):return False
    for annotation in chain['annotations']:
        if annotation.get('consumed_atoms'):
            if any(site_matches(a['site'],site,length) and wanted&set(a['atoms']) for a in annotation['consumed_atoms']):return False
        elif site in annotation.get('consumed_sites',[]):return False
    return not any(e['chain']==chain['id'] and site_matches(e['site'],site,length) and
                   wanted&atom_options(e['atom'],residue,e['site']) for b in state['bonds'] for e in b['endpoints'])

def endpoint_choices(state,chain,selector,atom,allowed_groups=()):
    p=chain['peptide']
    if 'nterm' in selector and p['nterm']!=selector['nterm'] or 'cterm' in selector and p['cterm']!=selector['cterm']:return
    for index in selected_sites(selector,p['sequence']):
        site=selector.get('site_label',index+1)
        if selector.get('group') and not any(m['site']==site and m.get('group')==selector['group'] for m in p['modifications']):continue
        if endpoint_available(state,chain,site,atom,allowed_groups):yield {'chain':chain['id'],'site':site,'atom':atom}

def unique_chain_id(state,base):
    ids={c['id'] for c in state['chains']};cid=base;i=1
    while cid in ids:cid=base+f'.{i}';i+=1
    return cid

def clone_chain(state,chain,cid):
    output=deepcopy(chain);output['id']=cid
    output['origin_sites']=[o.replace(chain['id']+':',cid+':',1) for o in output['origin_sites']]
    for location in output.get('charge_locations',[]):location['origin']=location['origin'].replace(chain['id']+':',cid+':',1)
    return output

def connection_transitions(state,record,branch,context):
    recipe=branch['recipe'];t=recipe['transform'];kind=t['kind'];selector=recipe['site_selector']
    allowed=t.get('allowed_occupied_groups',[])
    if kind=='trp_crosslink':source_atom='indole_C2';partners={'indole_C2':selector};source_site=None
    elif kind=='o_acylate':source_atom='sidechain_O';partners={};source_site=None
    elif kind=='n_acylate':source_atom=t['source_atom'];partners={};source_site=t.get('source_site')
    else:source_atom=t['source_atom'];partners=t['partner_atoms'];source_site=t.get('source_site')
    for source_chain in state['chains']:
        if not set(selector.get('requires_residues',[]))<=set(source_chain['peptide']['sequence']):continue
        source_selector={**selector,**({'site_label':source_site} if source_site is not None else {})}
        for left in endpoint_choices(state,source_chain,source_selector,source_atom,allowed):
            pool=[(state,c) for c in state['chains']]
            if kind in {'o_acylate','n_acylate'}:
                pool=[]
                for donor in context[0].get('acyl_donors',[]):
                    cid=unique_chain_id(state,source_chain['id']+'.acyl_donor')
                    c=normalize_molecule({'chains':[{'id':cid,**donor}],'bonds':[]})['chains'][0]
                    augmented={**state,'chains':state['chains']+[c]};pool.append((augmented,c))
                current_partners={'carbonyl_C':{'position':'end','site_label':'C','cterm':'acid'}}
            elif t.get('external_partner'):
                cid=unique_chain_id(state,source_chain['id']+'.external')
                c=normalize_molecule({'chains':[{'id':cid,**t['external_partner']}],'bonds':[]})['chains'][0]
                augmented={**state,'chains':state['chains']+[c]};pool=[(augmented,c)];current_partners=partners
            else:
                current_partners=partners
                # Homomolecular hypotheses use another explicit copy only when the source chain is an isolated molecule.
                if not any(e['chain']==source_chain['id'] for b in state['bonds'] for e in b['endpoints']):
                    cid=unique_chain_id(state,source_chain['id']+'.copy')
                    c=clone_chain(state,source_chain,cid)
                    pool.append(({**state,'chains':state['chains']+[c]},c))
            for augmented,partner in pool:
                for atom,partner_selector in current_partners.items():
                    for right in endpoint_choices(augmented,partner,partner_selector,atom,allowed):
                        intra=left['chain']==right['chain']
                        if left==right:continue
                        if intra and kind=='trp_crosslink' and not t.get('intrachain'):continue
                        if t.get('topologies') and ('intrachain' if intra else 'interchain') not in t['topologies']:continue
                        if t.get('symmetric') and encoded(left)>=encoded(right):continue
                        output=deepcopy(augmented)
                        delta={'H':-2,'O':-1} if kind in {'o_acylate','n_acylate'} else t['bond_delta']
                        bond={'kind':'O_acyl_ester' if kind=='o_acylate' else 'amide' if kind=='n_acylate' else t['bond'],
                              'endpoints':[deepcopy(left),deepcopy(right)],'delta_formula':deepcopy(delta)}
                        if kind=='trp_crosslink' and not t.get('symmetric'):bond['reduced_endpoint']=deepcopy(right)
                        if t.get('source_unsaturation'):bond['consumed_unsaturation']=deepcopy(left)
                        output['bonds'].append(bond);output['bonds'].sort(key=encoded);output['chains'].sort(key=lambda c:c['id'])
                        validate_state(output)
                        fit=conditions(recipe,local_context(context,source_chain))
                        yield {'rule_id':record['rule_id'],'branch_id':branch['branch_id'],
                            'sites':[site_origin(source_chain,left['site']),site_origin(partner,right['site'])],
                            'product_state':output,'delta_formula':add(state_formula(output),{e:-n for e,n in state_formula(state).items()}),
                            'contract':branch['contract'],'source':recipe['source_refs'],'reference_state':recipe['reference_state'],
                            'product':recipe['product_description'],'evidence_level':recipe['evidence_level'],**fit}

def network_transitions(state,record,branch,context):
    import itertools,math
    recipe=branch['recipe'];ss=[b for b in state['bonds'] if b['kind']=='disulfide'];count=len(ss)
    if not count:return
    fixed=[b for b in state['bonds'] if b['kind']!='disulfide']
    free_state={**state,'bonds':fixed};vertices=[]
    for chain in state['chains']:
        vertices.extend(endpoint_choices(free_state,chain,{'position':'any','residues':['C']},'Cys_Sgamma'))
    if len(vertices)<count*2:return
    total=math.comb(len(vertices),count*2)*math.factorial(count*2)//(2**count*math.factorial(count))
    if total>recipe['transform']['max_network_states']:raise UnsupportedTransition(f'Explicit SS network requires {total} states; the declared network budget is insufficient')
    def matchings(items):
        if not items:yield [];return
        for i,other in enumerate(items[1:],1):
            for rest in matchings(items[1:i]+items[i+1:]):yield [(items[0],other)]+rest
    original={encoded(sorted(b['endpoints'],key=encoded)) for b in ss}
    for chosen in itertools.combinations(vertices,count*2):
        for pairs in matchings(list(chosen)):
            if {encoded(sorted(pair,key=encoded)) for pair in pairs}==original:continue
            output=deepcopy(state)
            output['bonds']=deepcopy(fixed)+[{'kind':'disulfide','endpoints':sorted(deepcopy(pair),key=encoded),'delta_formula':{'H':-2}} for pair in pairs]
            output['bonds'].sort(key=encoded);validate_state(output)
            if state_formula(output)!=state_formula(state):raise ValueError('Conserved SS count changed whole-state composition')
            yield {'rule_id':record['rule_id'],'branch_id':branch['branch_id'],
                'sites':[e['chain']+':'+str(e['site']) for e in chosen],'product_state':output,'delta_formula':{},
                'event_kind':'network_endpoint_hypothesis','elementary_exchange_steps':'not_enumerated',
                'contract':branch['contract'],'source':recipe['source_refs'],'reference_state':recipe['reference_state'],
                'product':recipe['product_description'],'evidence_level':recipe['evidence_level'],**conditions(recipe,context)}
