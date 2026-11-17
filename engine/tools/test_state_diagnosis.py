import unittest
from pyteomics import mass
from state_diagnosis import diagnose_state_case
from test_structure_state import chain

def case(target,observations,rules,depth=2):
    return {'intended_product':target,'reaction':{'description':'state diagnostic test','stage_tags':['oxidation']},
            'observations':observations,'rule_ids':rules,'max_events':depth,'max_nodes':10000}
class StateDiagnosisTests(unittest.TestCase):
    def test_legacy_trp_and_cys_mass_cases_are_preserved(self):
        for seq,delta,rid in [('W',171.962827,'CH03-R009'),('GC',51.101428,'CH02-R002')]:
            c=case(chain(seq),[{'kind':'delta_mono','value':delta,'tolerance':0.00001,'tolerance_unit':'Da'}],[rid])
            c['reaction']['stage_tags']=[]
            r=diagnose_state_case(c)
            self.assertTrue(r['observations'][0]['candidates'])
            self.assertTrue(all(h['condition']=='unknown' for h in r['observations'][0]['candidates']))
    def test_pyroglutamate_closes_n_terminal_resource(self):
        c=case(chain('QG'),[{'kind':'delta_mono','value':-17.02654910112,'tolerance':0.00001,'tolerance_unit':'Da'}],['SUP-005','CH07-R013'])
        r=diagnose_state_case(c)
        self.assertTrue(r['observations'][0]['candidates'])
        self.assertFalse(any([e['rule_id'] for e in p]==['SUP-005','CH07-R013'] for h in r['observations'][0]['candidates'] for p in h['paths']))
    def test_endpoint_and_continuation_paths_share_one_candidate(self):
        c=case(chain('GC'),[{'kind':'delta_mono','value':31.98982923914,'reference':'intended_product','tolerance':1e-5,'tolerance_unit':'Da'}],
               ['CH09-R001','CH09-R002'])
        r=diagnose_state_case(c);hits=r['observations'][0]['candidates']
        self.assertEqual(len(hits),1)
        self.assertEqual({tuple(e['rule_id'] for e in p) for p in hits[0]['paths']},{('CH09-R002',),('CH09-R001','CH09-R002')})
    def test_two_real_peaks_with_different_charges(self):
        base=mass.calculate_mass(sequence='M');proton=1.007276466621
        obs={'kind':'peak_pair','reference':'intended_product','peaks':[
            {'kind':'mz','value':(base+2*proton)/2,'charge':2,'ions':{'H':2},'tolerance':0.001,'tolerance_unit':'Da'},
            {'kind':'mz','value':base+proton+15.99491461957,'charge':1,'ions':{'H':1},'tolerance':0.002,'tolerance_unit':'Da'}]}
        r=diagnose_state_case(case(chain('M'),[obs],['CH09-R029']))
        self.assertTrue(r['observations'][0]['candidates'])
        for h in r['observations'][0]['candidates']:
            self.assertAlmostEqual(h['propagated_error']['tolerance_da'],0.004,places=9)
            self.assertAlmostEqual(h['predicted'],15.99491461957,places=7)
    def test_declared_ion_hypotheses_and_unknown_conditions(self):
        observed=mass.calculate_mass(sequence='M')+15.99491461957+1.007276466621
        c=case(chain('M'),[{'kind':'mz','value':observed,'tolerance':0.00001,'tolerance_unit':'Da',
            'ion_hypotheses':[{'charge':1,'ions':{'H':1}},{'charge':1,'ions':{'Na':1}},{'charge':2,'ions':{'H':2}}]}],['CH09-R029'])
        c['reaction']['stage_tags']=[]
        r=diagnose_state_case(c);o=r['observations'][0]
        self.assertEqual(o['ion_hypothesis_count'],3);self.assertTrue(o['candidates'])
        self.assertTrue(all(h['condition']=='unknown' for h in o['candidates']))
    def test_explicit_conflicts_are_not_silently_dropped(self):
        c=case(chain('M'),[{'kind':'delta_mono','value':15.99491461957,'tolerance':0.00001,'tolerance_unit':'Da'}],['CH09-R029'])
        c['reaction']['stage_tags']=[]
        c['reaction']['stage_tags_absent']=['oxidation','air_oxidation','peroxide','h2o2','dmso','ozone','photooxidation']
        r=diagnose_state_case(c);o=r['observations'][0]
        self.assertFalse(o['candidates']);self.assertTrue(o['condition_conflicts'])
        self.assertTrue(any(e['conflicting_conditions'] for h in o['condition_conflicts'] for p in h['paths'] for e in p))
    def test_absolute_mass_after_interchain_reduction_uses_fragment_not_mass_sum(self):
        target={'chains':[chain('CG','a'),chain('GCA','b')],
                'bonds':[{'kind':'disulfide','endpoints':[{'chain':'a','site':1},{'chain':'b','site':2}]}]}
        obs={'kind':'neutral_mono','value':mass.calculate_mass(sequence='CG'),'tolerance':0.00001,'tolerance_unit':'Da'}
        r=diagnose_state_case(case(target,[obs],['CH09-R058'],1))
        h=r['observations'][0]['candidates'];self.assertTrue(h)
        self.assertTrue(all(x['component']['chain_ids']==['a'] for x in h))
    def test_heterogeneous_interchain_formation_keeps_chain_ids(self):
        target={'chains':[chain('CG','a'),chain('GCA','b')],'bonds':[]}
        value=mass.calculate_mass(sequence='CG')+mass.calculate_mass(sequence='GCA')-mass.calculate_mass(formula='H2')
        r=diagnose_state_case(case(target,[{'kind':'neutral_mono','value':value,'tolerance':1e-5,'tolerance_unit':'Da'}],['CH09-R005'],1))
        self.assertTrue(any(h['component']['chain_ids']==['a','b'] for h in r['observations'][0]['candidates']))
    def test_all_two_cys_endpoint_states_agree_with_independent_product_enumeration(self):
        # Independent oracle: each Cys is SH/SOH/SO2H/SO3H; no engine transform function used.
        base=mass.calculate_mass(sequence='CC');observations=[{'kind':'neutral_mono','value':base+n*15.99491461957,
                'tolerance':1e-5,'tolerance_unit':'Da'} for n in range(7)]
        r=diagnose_state_case(case(chain('CC'),observations,['CH09-R001','CH09-R002','CH09-R003']))
        lookup={'CysSulfenic':1,'CysSulfinic':2,'CysSulfonic':3}
        actual=set()
        for obs in r['observations']:
            for h in obs['candidates']:
                mods=h['component']['state']['chains'][0]['peptide']['modifications']
                values={m['site']:lookup[m['group']] for m in mods};actual.add((values.get(1,0),values.get(2,0)))
        expected={(a,b) for a in range(4) for b in range(4)}
        self.assertEqual(actual,expected)
        self.assertTrue(r['complete_for_compiled_transition_pool'])
