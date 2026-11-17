"""Independent, basic experimental helpers; predictions and measurements stay separate."""
from collections import Counter
import math
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'engine/tools'))
from case_input import read_sequence,check_keys
from structure_state import normalize_molecule,state_formula
from observation_model import formula_mass

def number(value,label,positive=False):
    if type(value) not in (int,float) or not math.isfinite(value) or (positive and value<=0):
        raise ValueError(label+' requires a finite '+('positive ' if positive else '')+'number')
    return value

def spps_review(raw):
    check_keys(raw,{'peptide','process','notes'},'SPPS review')
    peptide=raw.get('peptide',{});sequence=read_sequence(peptide.get('sequence'));process=raw.get('process',{})
    if not isinstance(process,dict):raise ValueError('process requires an object')
    check_keys(process,{'chemistry','resin','loading_mmol_g','coupling','deprotection','cleavage','conditions'},'SPPS process')
    if 'loading_mmol_g' in process:number(process['loading_mmol_g'],'resin loading',True)
    present=' '.join(str(v) for v in process.values()).lower();counts=Counter(sequence);findings=[]
    def flag(name,sites,reason,source='SPPS-REVIEW'):
        findings.append({'name':name,'sites':sites,'reason':reason,'source_id':source,'classification':'mechanistic_attention; not a failure probability'})
    for i,residue in enumerate(sequence):
        if residue=='D' and i+1<len(sequence):
            if sequence[i+1]=='P':flag('Asp-Pro酸解位点',[i+1,i+2],'核对酸接触及Asp-Pro断裂；不按普通Asp-X的下一肽键NH进攻模型归类。','YANG-2016')
            else:flag('Asp-X基序',[i+1,i+2],'核对碱暴露、保护态及aspartimide风险；不同X与工艺不能视为等概率。','BEHRENDT-2016')
        if residue in {'C','H'}:flag('消旋/侧链状态核对',[i+1],'活化、碱、温度及保护基会影响；这里只提示核对，不预测ee或dr。','YANG-2016')
        if residue=='M':flag('Met氧化敏感位点',[i+1],'核对氧化剂、空气/过氧化物及切割清除剂。','YANG-2016')
    if sequence[0] in {'Q','E'}:flag('N端环化参照',[1],'确认N端实际是否游离，再检查焦谷氨酸等环化；最终序列不等于树脂上的保护态。','YANG-2016')
    runs=[];start=0
    for i in range(1,len(sequence)+1):
        if i==len(sequence) or sequence[i]!=sequence[start]:
            if i-start>=3:runs.append({'residue':sequence[start],'start':start+1,'end':i,'length':i-start})
            start=i
    if runs:flag('重复片段/聚集核对',[r['start'] for r in runs],'重复片段是检查溶胀和聚集的线索；不把自由肽疏水性直接代入受保护肽-树脂体系。','SPPS-REVIEW')
    if 'methanol' in present or '甲醇' in present:
        if any(s in present for s in ('tfa','甲酸','酸')):flag('酸/甲醇处理',[i+1 for i,x in enumerate(sequence) if x in {'D','E'}],'检查侧链及C端甲酯；需按实际端基解释。','YANG-2016')
    mass=None
    missing=[x for x in ('nterm','cterm') if x not in peptide]
    if not missing:mass=formula_mass(state_formula(normalize_molecule(peptide)))
    return {'module':'spps_basics','sequence':sequence,'length':len(sequence),'residue_counts':dict(counts),
        'repeat_runs':runs,'reference_element_mass':mass,'missing_terminal_fields':missing,
        'missing_process_fields':[x for x in ('chemistry','resin','loading_mmol_g','coupling','deprotection','cleavage') if x not in process],
        'findings':findings,'process':process,'notes':raw.get('notes'),
        'scope':'sequence/process attention list; protected-residue behaviour, yield, purity and identity require separate experimental evidence',
        'predicted_success_probability':None}

def column_volume(column):
    check_keys(column,{'length_mm','id_mm','particle_um','pore_A','pH_range','name'},'HPLC column')
    length=number(column.get('length_mm'),'column length',True);diameter=number(column.get('id_mm'),'column ID',True)
    return math.pi*(diameter/2)**2*length/1000

def hplc_review(raw):
    check_keys(raw,{'purpose','peptide','column','flow_ml_min','gradient','mobile_phase','transfer_to','notes'},'HPLC method')
    column=raw.get('column',{});volume=column_volume(column);flow=number(raw.get('flow_ml_min'),'flow',True)
    phase=raw.get('mobile_phase',{});check_keys(phase,{'aqueous','organic','additives','pH'},'mobile phase')
    gradient=raw.get('gradient',[])
    if not isinstance(gradient,list) or len(gradient)<2:raise ValueError('Provide at least two gradient time/%B points')
    points=[]
    for point in gradient:
        check_keys(point,{'time_min','percent_B'},'gradient point');t=number(point.get('time_min'),'gradient time');b=number(point.get('percent_B'),'%B')
        if t<0 or not 0<=b<=100 or points and t<=points[-1]['time_min']:raise ValueError('Gradient times must increase and %B must be in 0..100')
        points.append({'time_min':t,'percent_B':b})
    slopes=[(b['percent_B']-a['percent_B'])/(b['time_min']-a['time_min']) for a,b in zip(points,points[1:])]
    additives=phase.get('additives',[])
    if not isinstance(additives,list) or any(not isinstance(x,str) for x in additives):raise ValueError('additives requires a string list')
    text=' '.join([str(phase.get('organic','')),*additives]).lower();findings=[]
    if 'tfa' in text or '三氟乙酸' in text:findings.append({'name':'TFA/ESI兼容性','reason':'LC分离需要与ESI响应分别评价；离子对添加剂可能抑制响应。','source_id':'HPLC-BASICS'})
    if any(x in text for x in ('phosphate','磷酸','sds')) and 'ms' in str(raw.get('purpose','')).lower():findings.append({'name':'非挥发添加剂','reason':'核对MS兼容性和仪器清洗要求。','source_id':'HPLC-METHOD'})
    if ('methanol' in text or '甲醇' in text or 'meoh' in text) and any(x in text for x in ('tfa','formic','acid','甲酸','酸')):
        findings.append({'name':'酸/甲醇甲酯化线索','reason':'检查D/E侧链及游离C端的处理副反应；不能把方法引入的峰当原始样品成分。','source_id':'YANG-2016'})
    if 'pH' in phase:
        pH=number(phase['pH'],'pH');limits=column.get('pH_range')
        if limits is not None:
            if not isinstance(limits,list) or len(limits)!=2 or limits[0]>limits[1]:raise ValueError('pH_range requires ordered lower/upper values')
            if not limits[0]<=pH<=limits[1]:findings.append({'name':'超出声明的柱pH范围','reason':'当前pH与输入的柱范围冲突。','source_id':'HPLC-METHOD'})
    transfer=None
    if raw.get('transfer_to'):
        new=raw['transfer_to'];newvolume=column_volume(new);newflow=flow*(new['id_mm']/column['id_mm'])**2
        scale=newvolume/volume*flow/newflow
        transfer={'suggested_flow_ml_min':newflow,'time_scale':scale,
            'gradient':[{'time_min':p['time_min']*scale,'percent_B':p['percent_B']} for p in points],
            'basis':'constant linear velocity and gradient column volumes; dwell volume, particle size, viscosity, pressure and selectivity still need verification'}
    return {'module':'hplc_basics','purpose':raw.get('purpose'),'geometric_column_volume_ml':volume,
        'flow_ml_min':flow,'gradient':points,'gradient_slopes_percent_B_per_min':slopes,'transfer':transfer,
        'findings':findings,'checklist':['实际柱化学/孔径/耐压','样品溶剂和进样负载','关键杂质对分辨率','MS添加剂兼容性','清洗/再平衡与空白','检测波长、纯度与回收率口径'],
        'scope':'review and scaling arithmetic; no validated retention or resolution prediction'}

def cyclic_review(raw):
    check_keys(raw,{'smiles','closure','expected_formula','expected_formal_charge','peptide_reference','reference_relation','notes'},'cyclic check')
    from rdkit import Chem,rdBase
    from rdkit.Chem import rdMolDescriptors
    with rdBase.BlockLogs():mol=Chem.MolFromSmiles(raw.get('smiles',''))
    if mol is None or mol.GetNumAtoms()==0:raise ValueError('A valid explicit SMILES is required')
    formula=rdMolDescriptors.CalcMolFormula(mol);charge=Chem.GetFormalCharge(mol);counts=Counter()
    for atom in mol.GetAtoms():counts[atom.GetSymbol()]+=1;counts['H']+=atom.GetTotalNumHs()
    counts={e:n for e,n in counts.items() if n};closure=raw.get('closure');check=None
    if closure:
        check_keys(closure,{'kind','atom_maps'},'closure');maps=closure.get('atom_maps')
        if not isinstance(maps,list) or len(maps)!=2 or any(type(x) is not int or x<=0 for x in maps) or len(set(maps))!=2:raise ValueError('Two distinct positive atom-map IDs are required')
        atoms={}
        for atom in mol.GetAtoms():
            n=atom.GetAtomMapNum()
            if n:
                if n in atoms:raise ValueError('Atom map IDs must be unique')
                atoms[n]=atom
        if not set(maps)<=set(atoms):raise ValueError('Declared closure maps missing from SMILES')
        a,b=[atoms[x] for x in maps];bond=mol.GetBondBetweenAtoms(a.GetIdx(),b.GetIdx());kind=closure.get('kind')
        if kind not in {'amide','disulfide'}:raise ValueError('Basic closure supports amide or disulfide')
        chemical=False
        if bond and bond.GetBondType()==Chem.BondType.SINGLE:
            if kind=='disulfide':chemical=a.GetSymbol()==b.GetSymbol()=='S'
            else:
                c=a if a.GetSymbol()=='C' else b if b.GetSymbol()=='C' else None;n=b if c is a else a
                chemical=c is not None and n.GetSymbol()=='N' and any(x.GetBondType()==Chem.BondType.DOUBLE and x.GetOtherAtom(c).GetSymbol()=='O' for x in c.GetBonds())
        check={'kind':kind,'atom_maps':maps,'bond_present':bond is not None,'bond_type_matches':chemical,
            'bond_in_ring':bool(bond and bond.IsInRing()),'pass':bool(chemical and bond.IsInRing())}
    formula_check=None
    if raw.get('expected_formula'):
        from mass_ledger import parse_formula
        formula_check=dict(parse_formula(raw['expected_formula']))==counts
    relation=None
    if raw.get('reference_relation'):
        if raw['reference_relation']!='head_to_tail':raise ValueError('Only the head_to_tail reference relation is implemented')
        ref=raw.get('peptide_reference',{})
        if ref.get('nterm')!='free' or ref.get('cterm')!='acid':raise ValueError('Head-to-tail water-loss comparison requires free N and acidic C termini')
        base=state_formula(normalize_molecule(ref));expected={e:base.get(e,0)-({'H':2,'O':1}.get(e,0)) for e in set(base)|{'H','O'}}
        expected={e:n for e,n in expected.items() if n};relation={'net_formula':{'H':-2,'O':-1},'matches_actual_formula':expected==counts and charge==0}
    expected_charge=raw.get('expected_formal_charge')
    if expected_charge is not None and type(expected_charge) is not int:raise ValueError('expected_formal_charge requires an integer')
    return {'module':'cyclic_basics','formula':formula,'element_counts':counts,'formal_charge':charge,
        'rdkit_exact_mass':rdMolDescriptors.CalcExactMolWt(mol),'ring_count':mol.GetRingInfo().NumRings(),
        'closure_check':check,'expected_formula_matches':formula_check,'expected_charge_matches':charge==expected_charge if expected_charge is not None else None,
        'reference_relation':relation,'scope':'2D graph/formula/formal-charge verification; permanent intrinsic charge, synthesis feasibility, conformational convergence and binding are not inferred'}
