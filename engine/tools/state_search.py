"""Enumerate bounded state paths. Completeness is false on truncation or unsupported transitions."""
from collections import deque
from copy import deepcopy
import time
from structure_state import encoded, fingerprint
from state_transformations import transitions

def path_key(path):
    return encoded(sorted([{'branch_id':e['branch_id'],'sites':e['sites'],
                            'reacting_endpoint':e.get('reacting_endpoint')} for e in path],key=encoded))

def binding_signature(state,event):
    resources=[]
    for source in event['sites']:
        found=False
        for chain in state['chains']:
            if source in chain['origin_sites']:
                i=chain['origin_sites'].index(source);p=chain['peptide'];found=True
                resources.append({'source':source,'residue':p['sequence'][i],
                    'neighbors':chain['origin_sites'][max(0,i-1):i+2],
                    'modifications':[m for m in p['modifications'] if m['site']==i+1],
                    'alpha':p.get('alpha_stereochemistry',{}).get(str(i+1)),
                    'relative_alpha':chain.get('relative_alpha_stereo',{}).get(source),
                    'bonds':[b for b in state['bonds'] if any(e['chain']==chain['id'] and e['site']==i+1 for e in b['endpoints'])]})
            elif source in {chain['id']+':N',chain['id']+':C'}:
                found=True;side=source.rsplit(':',1)[-1];p=chain['peptide']
                resources.append({'source':source,'terminus':p['nterm' if side=='N' else 'cterm'],
                                  'modifications':[m for m in p['modifications'] if m['site']==side],
                                  'terminal_origin':chain['origin_sites'][0 if side=='N' else -1]})
        if not found:resources.append({'source':source,'state':'explicit_partner_or_aggregate_binding'})
    return encoded({'branch_id':event['branch_id'],'resources':sorted(resources,key=encoded)})

def enumerate_states(bundle,initial,context,max_events=2,max_nodes=100000,rule_ids=None):
    if type(max_events) is not int or not 1<=max_events<=6:raise ValueError('Declare an event depth between 1 and 6')
    if type(max_nodes) is not int or max_nodes<1:raise ValueError('Node budget must be a positive integer')
    if rule_ids is not None and not set(rule_ids)<={r['rule_id'] for r in bundle['rules']}:raise ValueError('Unknown declared rule ID')
    start=time.perf_counter();nodes=[{'state':deepcopy(initial),'events':[]}]
    queue=deque([nodes[0]]);seen={(fingerprint(initial),'[]')};unsupported={};truncated=False;generated=0;transition_cache={}
    while queue:
        node=queue.popleft()
        if len(node['events'])>=max_events:continue
        cache_key=fingerprint(node['state'])
        if cache_key not in transition_cache:
            transition_cache[cache_key]=transitions(bundle,node['state'],context,None if rule_ids is None else set(rule_ids))
        events,gaps=transition_cache[cache_key]
        for gap in gaps:unsupported[encoded(gap)]=gap
        for event in events:
            generated+=1
            signature=binding_signature(node['state'],event)
            if any(e.get('input_binding_signature')==signature for e in node['events']):continue
            event={**event,'input_binding_signature':signature}
            path=node['events']+[event];state=event['product_state'];key=(fingerprint(state),path_key(path))
            if key in seen:continue
            if len(nodes)>=max_nodes:truncated=True;break
            seen.add(key);created={'state':state,'events':path};nodes.append(created);queue.append(created)
        if truncated:break
    selected=[r for r in bundle['rules'] if rule_ids is None or r['rule_id'] in rule_ids]
    no_templates=[r['rule_id'] for r in selected if not r['template']]
    unresolved={encoded(a):a for n in nodes for c in n['state']['chains'] for a in c['annotations'] if a.get('kind')=='unresolved_linkage'}
    return {'nodes':nodes,'scope':{'event_budget':max_events,'node_budget':max_nodes,'nodes':len(nodes),
        'events_generated':generated,'truncated':truncated,'truncation_reason':'node_budget' if truncated else None,
        'distinct_states_expanded':len(transition_cache),
        'unsupported_transitions':list(unsupported.values()),'source_records_without_templates':no_templates,
        'unresolved_product_structures':list(unresolved.values()),
        'complete_for_compiled_transition_pool':not truncated and not unsupported and not unresolved,
        'full_book_capabilities_complete':not truncated and not unsupported and not no_templates and not unresolved,
        'rule_ids':sorted(rule_ids) if rule_ids is not None else 'all_compiled_recipes',
        'seconds':time.perf_counter()-start}}
