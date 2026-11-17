import unittest
from knowledge_source import load_runtime
from structure_state import normalize_molecule,state_formula,components_of
from state_search import enumerate_states
from test_structure_state import chain

def context(*tags):return ({'stage_tags':set(tags),'features':set(),'residues':set()},
                           {'stage_tags':set(),'features':set(),'residues':set()})
class StateSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.bundle=load_runtime()[0]
    def search(self,state,rules,depth=2,tags=(),budget=10000):return enumerate_states(self.bundle,normalize_molecule(state),context(*tags),depth,budget,rules)
    def test_cys_oxidation_replaces_state_not_adds_endpoint_delta(self):
        r=self.search(chain('GC'),['CH09-R001','CH09-R002','CH09-R003'],3,['oxidation'])
        path=next(n for n in r['nodes'] if [e['rule_id'] for e in n['events']]==['CH09-R001','CH09-R002','CH09-R003'])
        self.assertEqual([e['delta_formula'] for e in path['events']],[{'O':1},{'O':1},{'O':1}])
        self.assertEqual(path['state']['chains'][0]['peptide']['modifications'],[{'site':2,'group':'CysSulfonic','name':'CysSulfonic'}])
        self.assertTrue(r['scope']['complete_for_compiled_transition_pool'])
    def test_deprotection_opens_a_site(self):
        r=self.search(chain('CG',modifications=[{'site':1,'group':'Acm'}]),['CH03-R021','CH09-R001'],2,['oxidation'])
        path=next(n for n in r['nodes'] if [e['rule_id'] for e in n['events']]==['CH03-R021','CH09-R001'])
        self.assertEqual(path['events'][1]['delta_formula'],{'O':1})
    def test_interchain_reduction_splits_components_and_opens_thiols(self):
        s={'chains':[chain('CG','a'),chain('GC','b')],'bonds':[{'kind':'disulfide','endpoints':[{'chain':'a','site':1},{'chain':'b','site':2}]}]}
        r=self.search(s,['CH09-R058','CH09-R001'],2,['oxidation'])
        path=next(n for n in r['nodes'] if [e['rule_id'] for e in n['events']]==['CH09-R058','CH09-R001'])
        self.assertEqual(len(components_of(path['state'])),2)
        self.assertEqual(path['events'][0]['delta_formula'],{'H':2})
    def test_two_positions_repeated_oxidation_are_not_collapsed(self):
        r=self.search(chain('CC'),['CH09-R001'],2,['oxidation'])
        doubles=[n for n in r['nodes'] if len(n['events'])==2]
        self.assertEqual(len(doubles),1)
        self.assertEqual({s for e in doubles[0]['events'] for s in e['sites']},{'chain1:1','chain1:2'})
    def test_budget_truncation_is_explicit(self):
        r=self.search(chain('CCC'),['CH09-R001','CH09-R002','CH09-R003'],3,['oxidation'],2)
        self.assertTrue(r['scope']['truncated']);self.assertFalse(r['scope']['complete_for_compiled_transition_pool'])
    def test_hydrolysis_keeps_both_coproduced_fragments(self):
        r=self.search(chain('ADPG'),['CH01-R008'],1)
        n=next(n for n in r['nodes'] if n['events'])
        self.assertEqual(sorted(''.join(c['peptide']['sequence']) for c in n['state']['chains']),['AD','PG'])
        self.assertEqual(n['events'][0]['delta_formula'],{'H':2,'O':1})
    def test_met_sulfoxide_to_sulfone_and_second_methyl_are_state_replacements(self):
        for sequence,rules,tags,expected in [('MG',['CH09-R029','CH09-R030'],['oxidation'],{'O':1}),
                                           ('KG',['CH07-R037','CH07-R038'],['formaldehyde','formic_acid'],{'C':1,'H':2})]:
            r=self.search(chain(sequence),rules,2,tags)
            paths=[n for n in r['nodes'] if [e['rule_id'] for e in n['events']]==rules and
                   n['events'][0]['sites']==n['events'][1]['sites']]
            self.assertTrue(paths)
            self.assertTrue(all(n['events'][1]['delta_formula']==expected for n in paths))
    def test_sequence_edit_remaps_existing_connections_by_residue_identity(self):
        s={'chains':[chain('HCG','a'),chain('GC','b')],
           'bonds':[{'kind':'disulfide','endpoints':[{'chain':'a','site':2},{'chain':'b','site':2}]}]}
        r=self.search(s,['CH05-R015'],1)
        n=next(n for n in r['nodes'] if n['events'] and n['events'][0]['sites']==['a:1'])
        a=next(c for c in n['state']['chains'] if c['id']=='a')
        self.assertEqual(a['peptide']['sequence'],['beta-Ala','H','C','G'])
        self.assertEqual(next(e['site'] for e in n['state']['bonds'][0]['endpoints'] if e['chain']=='a'),3)
