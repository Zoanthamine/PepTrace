"""Validate explicit neutral aldehyde/ketone components with the configured RDKit."""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from mass_ledger import MONO, parse_formula


def normalize_carbonyl_donors(items):
    if not isinstance(items, list):
        raise ValueError("carbonyl_donors requires an explicit component list")
    result = []
    for item in items:
        if not isinstance(item, dict) or set(item)-{"name","smiles"}:
            raise ValueError("A carbonyl donor permits only name and smiles")
        if not isinstance(item.get("smiles"), str) or not item["smiles"].strip():
            raise ValueError("A carbonyl donor requires explicit SMILES")
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol is None or len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol):
            raise ValueError("A carbonyl donor requires one valid neutral molecule")
        if any(a.GetIsotope() for a in mol.GetAtoms()):
            raise ValueError("Isotope-labelled carbonyl components require the future isotope model")
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol)))
        formula = parse_formula(rdMolDescriptors.CalcMolFormula(mol))
        if formula.keys()-MONO.keys():
            raise ValueError("Unsupported element in carbonyl component")
        ranks = list(Chem.CanonicalRankAtoms(mol,breakTies=False))
        sites = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            double_o = [b for b in atom.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(atom).GetAtomicNum()==8]
            if len(double_o) != 1:
                continue
            others = [n for n in atom.GetNeighbors() if n.GetIdx()!=double_o[0].GetOtherAtomIdx(atom.GetIdx())]
            if not all(n.GetAtomicNum()==6 for n in others) or atom.GetTotalNumHs()!=2-len(others):
                continue
            asymmetric = len(others)==1 or len(others)==2 and ranks[others[0].GetIdx()]!=ranks[others[1].GetIdx()]
            sites.append({"atom_index":atom.GetIdx(),"distinct_substituents":asymmetric})
        if not sites:
            raise ValueError("Carbonyl component must contain a ketone or aldehyde, not an acid/ester/amide")
        name = item.get("name",Chem.MolToSmiles(mol))
        if not isinstance(name,str) or not name.strip():
            raise ValueError("Carbonyl donor name must be nonempty")
        result.append({"name":name,"smiles":Chem.MolToSmiles(mol),"formula":formula,"carbonyl_sites":sites})
    return result
