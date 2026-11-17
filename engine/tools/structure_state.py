"""Canonical peptide chains, explicit connections, component masses and stable residue IDs."""
from copy import deepcopy
import hashlib
import json
from case_input import normalize_peptide, SULFUR_BOND_DELTAS
from knowledge_source import components
from mass_ledger import parse_formula, formula_mass
from mass_ledger import MONO
from atom_model import normalize_atom_peptide,canonical_atom,mod_atoms,atom_options,site_matches,atomic_occupied

CONNECTION_DELTAS={**SULFUR_BOND_DELTAS,'amide':{'H':-2,'O':-1},'O_acyl_ester':{'H':-2,'O':-1},
    'formaldehyde_methylene_bridge':{'C':1},'indole_indole':{'H':-2},'indole_indoline':{}}

def encoded(x):return json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(',',':'))
def add(*formulae):
    result={}
    for formula in formulae:
        for element,count in formula.items():result[element]=result.get(element,0)+count
    return {e:n for e,n in result.items() if n}

def chain_formula(peptide):
    c=components();table=c['residues']
    return add({'H':2,'O':1},*(table[r]['residue_formula'] for r in peptide['sequence']),
        {'C':2,'H':2,'O':1} if peptide['nterm']=='acetyl' else {},c['c_terminal_deltas'][peptide['cterm']],
        *(c['covalent_groups'][m['group']]['net_formula'] if 'group' in m else m['delta_formula'] for m in peptide['modifications']))

def normalize_molecule(raw):
    if 'chains' not in raw:
        raw={'chains':[{'id':'chain1',**deepcopy(raw)}],'bonds':[]}
    unknown=set(raw)-{'chains','bonds','geometry_status'}
    if unknown:raise ValueError('Unknown molecule fields: '+str(sorted(unknown)))
    if not isinstance(raw['chains'],list) or not raw['chains']:raise ValueError('Provide explicit peptide chains')
    chains=[];bonds=deepcopy(raw.get('bonds',[]));order=[]
    for i,entry in enumerate(raw['chains'],1):
        if not isinstance(entry,dict):raise ValueError('Chain requires an object')
        cid=entry.get('id',f'chain{i}')
        if not isinstance(cid,str) or not cid.strip() or cid in order:raise ValueError('Chain IDs must be unique nonempty strings')
        order.append(cid)
        charge=entry.get('intrinsic_charge',0)
        if type(charge) is not int:raise ValueError('Chain intrinsic charge must be an integer')
        peptide=normalize_atom_peptide({k:v for k,v in entry.items() if k not in {'id','intrinsic_charge','fixed_charge_sites'}})
        for bond in peptide.pop('bonds',[]):
            b={'kind':bond['kind'],'endpoints':[{'chain':cid,'site':site,'atom':'Cys_Sgamma'} for site in bond['sites']]}
            if 'oxidized_site' in bond:b['oxidized_endpoint']=next(e for e in b['endpoints'] if e['site']==bond['oxidized_site'])
            bonds.append(b)
        peptide.pop('topology',None)
        origins=[f'{cid}:{site}' for site in range(1,len(peptide['sequence'])+1)]
        locations=[]
        located_atoms=set()
        for location in entry.get('fixed_charge_sites',[]):
            if not isinstance(location,dict) or set(location)-{'site','atom','charge','name'}:raise ValueError('Unknown fixed-charge location fields')
            site=location.get('site');q=location.get('charge')
            if type(q) is not int or not q:raise ValueError('Fixed-charge location needs a nonzero integer charge')
            if type(site) is int:
                if not 1<=site<=len(origins):raise ValueError('Invalid fixed-charge site')
                origin=origins[site-1]
            elif site in {'N','C'}:origin=origins[0 if site=='N' else -1]
            else:raise ValueError('Invalid fixed-charge location site')
            atom=location.get('atom','Nalpha' if site=='N' else 'carbonyl_C' if site=='C' else None)
            if not isinstance(atom,str):raise ValueError('Fixed charge on a residue requires its explicit atom')
            atom=canonical_atom(atom,peptide['sequence'][origins.index(origin)],site)
            if (origin,atom) in located_atoms:raise ValueError('Repeated fixed-charge location at the same atom')
            located_atoms.add((origin,atom))
            locations.append({'origin':origin,'atom':atom,'charge':q,**({'name':location['name']} if 'name' in location else {})})
        localized=sum(x['charge'] for x in locations)
        if locations and 'intrinsic_charge' in entry and localized!=charge:raise ValueError('Charge locations do not add up to the declared intrinsic charge')
        if locations:charge=localized
        chain={'id':cid,'peptide':peptide,'origin_sites':origins,'intrinsic_charge':charge,'annotations':[]}
        if locations:chain['charge_locations']=locations
        elif charge:chain['unlocalized_intrinsic_charge']=charge
        chains.append(chain)
    byid={c['id']:c for c in chains};used=set();canonical=[]
    if not isinstance(bonds,list):raise ValueError('Bonds require a list')
    for bond in bonds:
        if set(bond)-{'kind','endpoints','delta_formula','oxidized_endpoint','lost_sulfur_endpoint'}:raise ValueError('Unknown connection fields')
        kind=bond['kind'];ends=deepcopy(bond['endpoints'])
        if not isinstance(ends,list) or len(ends)!=2:raise ValueError('A connection needs two endpoints')
        for e in ends:
            if set(e)-{'chain','site','atom'}:raise ValueError('Unknown endpoint fields')
            if type(e.get('chain')) is int:
                if not 1<=e['chain']<=len(order):raise ValueError('Unknown chain index')
                e['chain']=order[e['chain']-1]
            if e.get('chain') not in byid:raise ValueError('Unknown chain ID')
            chain=byid[e['chain']];site=e.get('site');peptide=chain['peptide']
            if type(site) is int:
                if not 1<=site<=len(peptide['sequence']):raise ValueError('Invalid connection site')
            elif site not in {'N','C'}:raise ValueError('Invalid connection site')
            atom=e.get('atom','Cys_Sgamma' if kind in SULFUR_BOND_DELTAS else None)
            if not isinstance(atom,str) or not atom:raise ValueError('An explicit connection atom is required')
            e['atom']=atom
            e['atom']=canonical_atom(atom,peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1],site)
            if kind in SULFUR_BOND_DELTAS and (type(site) is not int or peptide['sequence'][site-1]!='C'):
                raise ValueError('Sulfur bridges require two Cys sites')
            wanted=atom_options(e['atom'],peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1],site)
            if any(site_matches(m['site'],site,len(peptide['sequence'])) and wanted&mod_atoms(peptide,m) for m in peptide['modifications']):raise ValueError('Connection overlaps a modified atom')
            if site=='N' and peptide['nterm']!='free' or site=='C' and peptide['cterm']!='acid':raise ValueError('Connection overlaps a terminal state')
            resource=(e['chain'],site,atom)
            if resource in used:raise ValueError('Connection endpoint is consumed by another bond')
            used.add(resource)
        if ends[0]==ends[1]:raise ValueError('A connection cannot join an endpoint to itself')
        delta=bond.get('delta_formula',CONNECTION_DELTAS.get(kind))
        if delta is None:raise ValueError('Unsupported connection without explicit element bookkeeping: '+kind)
        if kind in CONNECTION_DELTAS and delta!=CONNECTION_DELTAS[kind]:raise ValueError('Connection formula contradicts the defined transformation')
        if not isinstance(delta,dict) or any(e not in MONO or type(n) is not int for e,n in delta.items()):
            raise ValueError('Invalid connection elemental delta')
        output={'kind':kind,'endpoints':sorted(ends,key=encoded),'delta_formula':delta}
        if kind in {'thiosulfinate','thiosulfonate'}:
            ox=bond.get('oxidized_endpoint')
            if ox is None:raise ValueError('Oxidized bridge needs an explicit oxidized endpoint')
            ox=deepcopy(ox)
            if type(ox.get('chain')) is int:ox['chain']=order[ox['chain']-1]
            ox.setdefault('atom','Cys_Sgamma')
            if ox not in ends:raise ValueError('Oxidized endpoint must belong to its bridge')
            output['oxidized_endpoint']=ox
        canonical.append(output)
    state={'chains':sorted(chains,key=lambda c:c['id']),'bonds':sorted(canonical,key=encoded),'geometry_status':'not_assessed'}
    validate_state(state)
    return state

def validate_state(state):
    ids=[c['id'] for c in state['chains']]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate product chain IDs')
    for chain in state['chains']:
        if len(chain['origin_sites'])!=len(chain['peptide']['sequence']):raise ValueError('Residue identity map is inconsistent')
        if any(n<0 for n in chain_formula(chain['peptide']).values()):raise ValueError('Negative product composition')
        localized=sum(x['charge'] for x in chain.get('charge_locations',[]))
        if localized+chain.get('unlocalized_intrinsic_charge',0)!=chain['intrinsic_charge']:
            raise ValueError('Product charge locations and intrinsic charge disagree')
        if any(x['origin'] not in chain['origin_sites'] for x in chain.get('charge_locations',[])):
            raise ValueError('Charge location lost its residue identity')
    byid={c['id']:c for c in state['chains']};used=set()
    for bond in state['bonds']:
        for endpoint in bond['endpoints']:
            if endpoint['chain'] not in byid:raise ValueError('Connection refers to a missing product chain')
            chain=byid[endpoint['chain']];site=endpoint['site'];p=chain['peptide']
            if type(site) is int and not 1<=site<=len(p['sequence']):raise ValueError('Connection site lost its residue identity')
            if type(site) is not int and site not in {'N','C'}:raise ValueError('Invalid product endpoint site')
            resource=(endpoint['chain'],site,endpoint['atom'])
            if resource in used:raise ValueError('A product atom was consumed by two connections')
            used.add(resource)
    if any(n<0 for n in state_formula(state).values()):raise ValueError('Negative connected molecular composition')

def state_formula(state):return add(*(chain_formula(c['peptide']) for c in state['chains']),*(b['delta_formula'] for b in state['bonds']))
def state_charge(state):return sum(c['intrinsic_charge'] for c in state['chains'])
def fingerprint(state):
    canonical=deepcopy(state)
    canonical['chains'].sort(key=lambda c:c['id']);canonical['bonds'].sort(key=encoded)
    for c in canonical['chains']:
        c['annotations'].sort(key=encoded)
        c['peptide']['modifications'].sort(key=encoded)
    return hashlib.sha256(encoded(canonical).encode()).hexdigest()

def physical_fingerprint(state):
    canonical=deepcopy(state)
    for chain in canonical['chains']:
        chain['annotations']=[a for a in chain['annotations'] if a.get('kind')!='product_provenance']
        for m in chain['peptide']['modifications']:
            if 'group' in m:m.pop('name',None)
    return fingerprint(canonical)

def components_of(state):
    left={c['id'] for c in state['chains']};groups=[]
    while left:
        chosen={min(left)}
        while True:
            updated=chosen|{e['chain'] for b in state['bonds'] if chosen&{e['chain'] for e in b['endpoints']} for e in b['endpoints']}
            if updated==chosen:break
            chosen=updated
        left-=chosen
        component={'chains':[deepcopy(c) for c in state['chains'] if c['id'] in chosen],
                   'bonds':[deepcopy(b) for b in state['bonds'] if {e['chain'] for e in b['endpoints']}<=chosen],
                   'geometry_status':'not_assessed'}
        groups.append(component)
    return groups

def projection_formula(state,cid,peptide,recipe=None):
    chain=next(c for c in state['chains'] if c['id']==cid)
    if recipe is not None:return peptide['sequence'],chain_formula(peptide),atomic_occupied(state,chain,peptide,recipe)
    occupied={m['site'] for m in peptide['modifications']}
    occupied|={e['site'] for b in state['bonds'] for e in b['endpoints'] if e['chain']==cid}
    occupied|={site for annotation in chain['annotations'] for site in annotation.get('consumed_sites',[])}
    return peptide['sequence'],chain_formula(peptide),occupied

def describe(state):
    result=[]
    for component in components_of(state):
        formula=state_formula(component);q=state_charge(component)
        result.append({'state':component,'formula':formula,'neutral_mono':formula_mass(formula) if q==0 else None,
                       'element_mass':formula_mass(formula),'intrinsic_charge':q,'chain_ids':[c['id'] for c in component['chains']]})
    return result
