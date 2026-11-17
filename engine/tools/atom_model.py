"""Named reaction atoms and functional-slot occupancy; ambiguous atom assignments stay explicit."""
from copy import deepcopy
from case_input import normalize_peptide

BACKBONE={'Nalpha','backbone_N','alpha_C','beta_C','carbonyl_C','carbonyl_O'}
SIDE={
 'C':{'Cys_Sgamma'},'M':{'Met_Sdelta'},'K':{'Lys_Nepsilon'},'Orn':{'Orn_Ndelta'},
 'N':{'sidechain_amide_N','sidechain_carbonyl_C','sidechain_carbonyl_O'},
 'Q':{'sidechain_amide_N','sidechain_carbonyl_C','sidechain_carbonyl_O'},
 'D':{'sidechain_carboxyl_O','sidechain_carbonyl_C','sidechain_carbonyl_O'},
 'E':{'sidechain_carboxyl_O','sidechain_carbonyl_C','sidechain_carbonyl_O'},
 'S':{'sidechain_O'},'T':{'sidechain_O'},'Hyp':{'sidechain_O'},'Dap':{'sidechain_amine_N'},
 'H':{'imidazole_Npi','imidazole_Ntau','imidazole_C2','imidazole_C4','imidazole_C5'},
 'W':{'indole_N1','indole_C2','indole_C3','indole_C4','indole_C5','indole_C6','indole_C7'},
 'Y':{'phenolic_O','aryl_C2','aryl_C3','aryl_C5','aryl_C6'},
 'F':{'aryl_C2','aryl_C3','aryl_C4','aryl_C5','aryl_C6'},
 'Bpa':{'aryl_C2','aryl_C3','aryl_C4','aryl_C5','aryl_C6','ketone_C','ketone_O'},
 'R':{'guanidino_Ndelta','guanidino_Nomega','guanidino_Nomega_prime'},
 'Aoa':{'aminooxy_N','aminooxy_O'},'Dha':{'beta_C'},'MeDha':{'beta_C'}}
ALIASES={'indole_N':'indole_N1','phenol_ortho_C3':'aryl_C3','phenol_ortho_C5':'aryl_C5',
         'Gly_Nalpha':'Nalpha','side_chain_O':'sidechain_O','sidechain_N':'sidechain_amine_N','Dha_Cbeta':'beta_C'}
AMBIGUOUS={'site_local','imidazole_N','guanidino_N','indole_C','aryl_C'}

def canonical_atom(atom,residue,site):
    atom=ALIASES.get(atom,atom)
    if atom=='sidechain_O' and residue=='Y':atom='phenolic_O'
    if atom=='sidechain_O' and residue in {'D','E'}:atom='sidechain_carboxyl_O'
    if atom=='carbonyl_C' and site=='C':return atom
    if atom=='Nalpha' and site=='N':return atom
    if atom=='Nalpha' and site not in {'N',1}:raise ValueError('Internal peptide nitrogen must be identified as backbone_N')
    if atom=='backbone_N' and site in {'N',1}:return 'Nalpha'
    if atom=='cterm_hydrazide_Nprime' and site=='C':return atom
    allowed=BACKBONE|SIDE.get(residue,set())|AMBIGUOUS
    if residue in {'S','T','Y','Hyp'}:allowed=allowed|{'phosphate_P'}
    if atom not in allowed:raise ValueError(f'Atom {atom} does not belong to residue {residue}')
    return atom

def atom_options(atom,residue,site):
    atom=canonical_atom(atom,residue,site)
    if atom=='site_local':return SIDE.get(residue,{'site_local'})
    if atom=='imidazole_N':return {'imidazole_Npi','imidazole_Ntau'}
    if atom=='guanidino_N':return SIDE['R']
    if atom=='indole_C':return {a for a in SIDE['W'] if a.startswith('indole_C')}
    if atom=='aryl_C':return {a for a in SIDE.get(residue,set()) if a.startswith('aryl_C')}
    return {atom}

def default_mod_atom(peptide,mod):
    site=mod['site'];group=mod.get('group');name=mod.get('name','')
    if site=='N':return 'Nalpha'
    if site=='C':return 'cterm_hydrazide_Nprime' if peptide['cterm']=='hydrazide' and group not in {'C-amide','C-methyl-ester'} else 'carbonyl_C'
    residue=peptide['sequence'][site-1]
    if group=='Methyl' and residue=='M':raise ValueError('Methyl = CH2 cannot describe Met S-methyl sulfonium; declare CH3 and its fixed charge')
    if residue=='C':return 'Cys_Sgamma'
    if residue=='M':return 'Met_Sdelta'
    if residue in {'K','Orn','Dap'}:return {'K':'Lys_Nepsilon','Orn':'Orn_Ndelta','Dap':'sidechain_amine_N'}[residue]
    if residue in {'D','E'}:return 'sidechain_carboxyl_O'
    if residue in {'S','T','Hyp'}:return 'sidechain_O'
    if residue in {'N','Q'}:return 'sidechain_amide_N'
    if residue=='R':return 'guanidino_N'
    if residue=='H':return 'imidazole_N'
    if residue=='W':
        if group in {'Boc','Formyl','Cbz','Alloc'} or 'N1' in name or 'indole N' in name:return 'indole_N1'
        for number in range(2,8):
            if f'C{number}' in name:return f'indole_C{number}'
        return 'indole_C'
    if residue=='Y':return 'phenolic_O' if group in {'Ac','Formyl','tBu','All','Bzl','Phospho','Sulfo','Trifluoroacetyl'} else 'aryl_C'
    if residue=='F':return 'aryl_C'
    return 'site_local'

def mod_atoms(peptide,mod):
    site=mod['site'];residue=peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1]
    values=mod.get('atoms',[mod.get('atom') or default_mod_atom(peptide,mod)])
    if not isinstance(values,list) or not values or any(not isinstance(x,str) for x in values):raise ValueError('Modification atoms must be named')
    return {a for value in values for a in atom_options(value,residue,site)}

def site_matches(a,b,length):
    return a==b or (a=='N' and b==1 or b=='N' and a==1) or (a=='C' and b==length or b=='C' and a==length)

def normalize_atom_peptide(raw):
    """Validate separately located covalent groups instead of rejecting an entire residue."""
    mods=raw.get('modifications',[])
    if not isinstance(mods,list):raise ValueError('Modifications require a list')
    base={k:deepcopy(v) for k,v in raw.items() if k not in {'modifications','monomer_state'}}
    peptide=normalize_peptide({**base,'modifications':[]})
    result=[];occupied=[];monomer=None
    for original in mods:
        if not isinstance(original,dict):raise ValueError('Modification requires an object')
        mod=deepcopy(original);atom=mod.get('atom');atoms=mod.get('atoms')
        if atom is not None and atoms is not None:raise ValueError('Specify atom or atoms, not both')
        site=mod.get('site');length=len(peptide['sequence'])
        if type(site) is int:
            if not 1<=site<=length:raise ValueError('Invalid modification site')
            if atom=='Nalpha' or atom=='backbone_N' and site==1:
                if site!=1:raise ValueError('Internal backbone nitrogen must be called backbone_N')
                mod['site']='N';site='N'
            elif atom=='carbonyl_C' and site==length and mod.get('group') in {'C-amide','C-methyl-ester'}:mod['site']='C';site='C'
        core={k:v for k,v in mod.items() if k not in {'atom','atoms'}}
        single=normalize_peptide({**base,'modifications':[core],**({'monomer_state':raw['monomer_state']} if raw.get('monomer_state') else {})})
        if single.get('monomer_state'):monomer=single['monomer_state']
        resource_atoms=mod_atoms(peptide,mod)
        for previous,p_atoms in occupied:
            if site_matches(site,previous['site'],length) and resource_atoms&p_atoms:
                raise ValueError('Two modifications consume the same atom or unresolved functional slot')
            nh_family={'imidazole_Npi','imidazole_Ntau'}
            if (site_matches(site,previous['site'],length) and resource_atoms&nh_family and p_atoms&nh_family
                    and mod.get('group') and previous.get('group')):
                raise ValueError('Neutral His imidazole N-H substitutions share one available N-H; a second quaternizing step needs an explicit charge/formula')
        occupied.append((mod,resource_atoms))
        if single['nterm']!=base['nterm']:
            if peptide['nterm']!=base['nterm']:raise ValueError('Repeated N-terminal state')
            peptide['nterm']=single['nterm']
        elif single['cterm']!=base['cterm']:
            if peptide['cterm']!=base['cterm']:raise ValueError('Repeated C-terminal state')
            peptide['cterm']=single['cterm']
        else:
            core=single['modifications'][0]
            if atom is not None:core['atom']=canonical_atom(atom,peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1],site)
            if atoms is not None:core['atoms']=[canonical_atom(a,peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1],site) for a in atoms]
            result.append(core)
    peptide['modifications']=result
    if monomer:peptide['monomer_state']=monomer
    return peptide

def query_atoms(recipe,peptide,site):
    t=recipe['transform'];kind=t['kind'];selector=recipe.get('site_selector',{})
    if kind=='epimerize_alpha':return {'alpha_C'}
    if site=='N':return {'Nalpha'}
    if kind in {'n_terminal_delta','capped_suffix','carbonyl_condensation'}:
        if type(site) is int and site>1 and selector.get('next_backbone_amide'):return {'backbone_N'}
        if selector.get('side_chain_required') and type(site) is int:
            return {'Nalpha'}|mod_atoms(peptide,{'site':site,'name':t.get('modification','')})
        return {'Nalpha'}
    if kind in {'c_terminal_delta','c_terminal_acid','c_terminal_amide_substitution','c_terminal_pyrrolidinamide'} or site=='C':
        return {'cterm_hydrazide_Nprime' if peptide['cterm']=='hydrazide' and kind=='c_terminal_delta' else 'carbonyl_C'}
    if kind in {'delete','duplicate_residue','insert_before','insert_gly_before','prefix_dkp','hydrolytic_split','acyl_donor_self_condensation'}:
        return BACKBONE
    if kind=='backbone_amide_to_ester':return {'backbone_N','sidechain_O','carbonyl_C'}
    if selector.get('group'):
        relevant=[m for m in peptide['modifications'] if m['site']==site and m.get('group')==selector['group']]
        if relevant:return set.union(*(mod_atoms(peptide,m) for m in relevant))
    pseudo={'site':site,'group':t.get('product_group'),'name':t.get('modification','')}
    if t.get('reaction_atoms'):return set(t['reaction_atoms'])
    return mod_atoms(peptide,pseudo)

def modification_matches_selector(recipe,peptide,mod,site):
    selector=recipe.get('site_selector',{});wanted=selector.get('atom')
    if wanted:
        residue=peptide['sequence'][0 if site=='N' else -1 if site=='C' else site-1]
        return bool(atom_options(wanted,residue,site)&mod_atoms(peptide,mod))
    implied={'site':site,'group':mod.get('group')}
    return bool(mod_atoms(peptide,implied)&mod_atoms(peptide,mod))

def consumed_atoms(recipe,peptide,event):
    t=recipe['transform'];site=event['site'];cycle=t.get('local_cycle','')
    if type(site) is not int:site=next((x for x in event['consumed'] if type(x) is int),1)
    if cycle in {'oxazoline','oxazole'}:
        return [{'site':site-1,'atoms':['carbonyl_C']},{'site':site,'atoms':['backbone_N','sidechain_O']}]
    if cycle in {'aspartimide','glutarimide'} or '酰亚胺' in t.get('modification',''):
        return [{'site':site,'atoms':['sidechain_carbonyl_C','sidechain_carboxyl_O']},
                {'site':site+1,'atoms':['backbone_N']}]
    if t.get('consume_nterm') or cycle=='pyroglutamyl':
        return [{'site':'N','atoms':['Nalpha']},{'site':site,'atoms':['sidechain_carbonyl_C']}]
    result=[]
    for s in event['consumed']:
        if type(s) is int or s in ('N','C'):result.append({'site':s,'atoms':sorted(query_atoms(recipe,peptide,s))})
    return result

def atomic_occupied(state,chain,peptide,recipe):
    occupied=set();length=len(peptide['sequence'])
    for site in list(range(1,length+1))+['N','C']:
        wanted=query_atoms(recipe,peptide,site)
        if site=='C' and any(m.get('group') in {'NCA','NChloroformylAcidChloride'} for m in peptide['modifications']):occupied.add(site)
        for mod in peptide['modifications']:
            if site_matches(site,mod['site'],length) and wanted&mod_atoms(peptide,mod):occupied.add(site)
        for bond in state['bonds']:
            for endpoint in bond['endpoints']:
                if endpoint['chain']!=chain['id'] or not site_matches(site,endpoint['site'],length):continue
                residue=peptide['sequence'][0 if endpoint['site']=='N' else -1 if endpoint['site']=='C' else endpoint['site']-1]
                if wanted&atom_options(endpoint['atom'],residue,endpoint['site']):occupied.add(site)
        for annotation in chain['annotations']:
            explicit=annotation.get('consumed_atoms')
            if explicit:
                if any(site_matches(site,a['site'],length) and wanted&set(a['atoms']) for a in explicit):occupied.add(site)
            elif site in annotation.get('consumed_sites',[]):occupied.add(site)
        if type(site) is int and recipe.get('site_selector',{}).get('unmodified_side_chain'):
            side=SIDE.get(peptide['sequence'][site-1],set())
            if any(m['site']==site and side&mod_atoms(peptide,m) for m in peptide['modifications']):occupied.add(site)
    return occupied
