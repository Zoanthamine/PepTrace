import unittest
from pyteomics import mass
from structure_state import normalize_molecule,state_formula,components_of,state_charge

def chain(seq,cid=None,**extra):return {'sequence':seq,'nterm':'free','cterm':'acid',**({'id':cid} if cid else {}),**extra}
class StructureStateTests(unittest.TestCase):
    def test_heterogeneous_disulfide_and_chain_order(self):
        bond={'kind':'disulfide','endpoints':[{'chain':'a','site':1},{'chain':'b','site':2}]}
        a=normalize_molecule({'chains':[chain('CG','a'),chain('GC','b')],'bonds':[bond]})
        b=normalize_molecule({'chains':[chain('GC','b'),chain('CG','a')],'bonds':[bond]})
        self.assertEqual(state_formula(a),state_formula(b));self.assertEqual(len(components_of(a)),1)
        self.assertAlmostEqual(mass.calculate_mass(composition=state_formula(a)),
                               mass.calculate_mass(sequence='CG')+mass.calculate_mass(sequence='GC')-mass.calculate_mass(formula='H2'),places=8)
    def test_single_chain_bridge_equivalence(self):
        a=normalize_molecule(chain('CGC',bonds=[{'kind':'disulfide','sites':[1,3]}]))
        b=normalize_molecule({'chains':[chain('CGC','chain1')],'bonds':[{'kind':'disulfide','endpoints':[{'chain':1,'site':1},{'chain':1,'site':3}]}]})
        self.assertEqual(a,b)
    def test_illegal_occupancy_and_formula_are_rejected(self):
        for raw in [
            {'chains':[chain('CG','a'),chain('GC','b')],'bonds':[{'kind':'disulfide','endpoints':[{'chain':'a','site':1},{'chain':'b','site':1}]}]},
            {'chains':[chain('CC','a',modifications=[{'site':1,'group':'Acm'}])],'bonds':[{'kind':'disulfide','endpoints':[{'chain':'a','site':1},{'chain':'a','site':2}]}]},
            {'chains':[chain('CC','a')],'bonds':[{'kind':'disulfide','delta_formula':{'H':-1},'endpoints':[{'chain':'a','site':1},{'chain':'a','site':2}]}]},
        ]:
            with self.subTest(raw=raw),self.assertRaises(ValueError):normalize_molecule(raw)
    def test_disconnected_charges_are_kept_per_component(self):
        state=normalize_molecule({'chains':[chain('G','a',intrinsic_charge=1),chain('A','b')],'bonds':[]})
        self.assertEqual(sorted(state_charge(x) for x in components_of(state)),[0,1])
