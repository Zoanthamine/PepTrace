import unittest
from pyteomics import mass
from knowledge_source import load_runtime
from structure_state import normalize_molecule,chain_formula,components_of,state_charge
from state_search import enumerate_states
from state_diagnosis import diagnose_state_case
from test_structure_state import chain
from test_state_search import context

class AtomChargeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.bundle=load_runtime()[0]
    def search(self,raw,rules,depth=1):return enumerate_states(self.bundle,normalize_molecule(raw),context(),depth,10000,rules)
    def test_sidechain_reaction_remains_possible_with_n_terminal_protection(self):
        raw=chain('KG',modifications=[{'site':'N','group':'Boc'}])
        result=self.search(raw,['CH07-R003'])
        product=next(n for n in result['nodes'] if n['events'])
        self.assertEqual(product['events'][0]['delta_formula'],{'C':2,'H':2,'O':1})
        self.assertTrue(any(m.get('group')=='Boc' for m in product['state']['chains'][0]['peptide']['modifications']))
    def test_trp_n_boc_does_not_occupy_c2(self):
        result=self.search(chain('WG',modifications=[{'site':1,'group':'Boc','atom':'indole_N1'}]),['CH03-R001'])
        self.assertTrue(any(n['events'] for n in result['nodes']))
    def test_two_distinct_atoms_on_one_residue_are_allowed(self):
        raw=chain('ACKG',modifications=[{'site':2,'group':'Methyl','atom':'backbone_N'},
                                       {'site':2,'group':'Acm','atom':'Cys_Sgamma'}])
        state=normalize_molecule(raw)
        self.assertEqual(len(state['chains'][0]['peptide']['modifications']),2)
        expected=mass.calculate_mass(sequence='ACKG')+mass.calculate_mass(formula='CH2')+mass.calculate_mass(formula='C3H5NO')
        self.assertAlmostEqual(mass.calculate_mass(composition=chain_formula(state['chains'][0]['peptide'])),expected,places=7)
        result=enumerate_states(self.bundle,state,context(),1,10000,['CH03-R021'])
        product=next(n for n in result['nodes'] if n['events'])
        remaining=product['state']['chains'][0]['peptide']['modifications']
        self.assertEqual(remaining,[{'site':2,'group':'Methyl','atom':'backbone_N'}])
    def test_duplicate_atom_and_incompatible_atom_are_rejected(self):
        for raw in [chain('K',modifications=[{'site':'N','group':'Boc'},{'site':1,'atom':'Nalpha','group':'Ac'}]),
                    chain('S',modifications=[{'site':1,'atom':'Cys_Sgamma','delta_formula':{'O':1}}]),
                    chain('M',modifications=[{'site':1,'group':'Methyl'}])]:
            with self.subTest(raw=raw),self.assertRaises(ValueError):normalize_molecule(raw)
    def test_his_neutral_n_substitution_capacity_is_not_duplicated(self):
        with self.assertRaises(ValueError):normalize_molecule(chain('H',modifications=[
            {'site':1,'group':'Formyl','atom':'imidazole_Npi'},
            {'site':1,'group':'Formyl','atom':'imidazole_Ntau'}]))
    def test_located_n_terminal_charge_stays_with_n_fragment(self):
        raw=chain('GADPG',intrinsic_charge=1,
            modifications=[{'site':'N','atom':'Nalpha','delta_formula':{'C':3,'H':7},'name':'trimethyl N'}],
            fixed_charge_sites=[{'site':'N','atom':'Nalpha','charge':1}])
        result=self.search(raw,['CH01-R008'])
        product=next(n for n in result['nodes'] if n['events'])
        components=components_of(product['state'])
        charges={''.join(c['chains'][0]['peptide']['sequence']):state_charge(c) for c in components}
        self.assertEqual(charges,{'GAD':1,'PG':0})
        expected=mass.calculate_mass(sequence='GAD')+mass.calculate_mass(formula='C3H7')-0.000548579909
        case={'intended_product':raw,'reaction':{'description':'explicit charged fragmentation'},'rule_ids':['CH01-R008'],'max_events':1,
              'observations':[{'kind':'mz','value':expected,'charge':1,'intrinsic_charge':1,'ions':{},'tolerance':1e-5,'tolerance_unit':'Da'}]}
        self.assertTrue(diagnose_state_case(case)['observations'][0]['candidates'])
    def test_located_sidechain_charge_stays_with_c_fragment(self):
        raw=chain('GADPK',intrinsic_charge=1,
            modifications=[{'site':5,'atom':'Lys_Nepsilon','delta_formula':{'C':3,'H':7}}],
            fixed_charge_sites=[{'site':5,'atom':'Lys_Nepsilon','charge':1}])
        result=self.search(raw,['CH01-R008']);product=next(n for n in result['nodes'] if n['events'])
        self.assertEqual({''.join(c['chains'][0]['peptide']['sequence']):state_charge(c) for c in components_of(product['state'])},{'GAD':0,'PK':1})
    def test_unlocated_charge_does_not_get_guessed_on_fragmentation(self):
        result=self.search(chain('GADPG',intrinsic_charge=1),['CH01-R008'])
        self.assertFalse(any(n['events'] for n in result['nodes']))
        self.assertTrue(result['scope']['unsupported_transitions'])
    def test_charge_sites_must_sum_to_declared_charge(self):
        with self.assertRaises(ValueError):normalize_molecule(chain('KG',intrinsic_charge=2,fixed_charge_sites=[{'site':1,'atom':'Lys_Nepsilon','charge':1}]))
    def test_charge_location_is_reindexed_after_insertion(self):
        raw=chain('HK',intrinsic_charge=1,fixed_charge_sites=[{'site':2,'atom':'Lys_Nepsilon','charge':1}])
        result=self.search(raw,['CH05-R015']);product=next(n for n in result['nodes'] if n['events'] and n['events'][0]['sites']==['chain1:1'])
        c=product['state']['chains'][0]
        self.assertEqual(c['origin_sites'].index(c['charge_locations'][0]['origin']),2)
        self.assertEqual(c['intrinsic_charge'],1)
