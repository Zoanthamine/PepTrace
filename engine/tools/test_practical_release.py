import unittest
from state_diagnosis import diagnose_state_case

class PracticalReleaseTests(unittest.TestCase):
    def test_unimplemented_structure_is_a_source_lead_not_an_executed_candidate(self):
        c={'intended_product':{'sequence':'GAGA','nterm':'free','cterm':'acid'},'reaction':{'description':'source-only family'},
           'observations':[{'kind':'delta_mono','value':41.02654910112,'tolerance':0.00001,'tolerance_unit':'Da'}],
           'max_events':1}
        r=diagnose_state_case(c);o=r['observations'][0]
        lead=next(x for x in o['source_leads'] if x['id']=='AI-054')
        self.assertFalse(lead['structure_confirmed']);self.assertEqual(lead['classification'],'source_lead_only')
        self.assertFalse(o['candidates'])
    def test_low_resolution_book_integer_is_not_a_high_resolution_confirmation(self):
        c={'intended_product':{'sequence':'GAGA','nterm':'free','cterm':'acid'},'reaction':{'description':'precision'},
           'observations':[{'kind':'delta_mono','value':41.5,'tolerance':0.00001,'tolerance_unit':'Da'}],
           'max_events':1}
        r=diagnose_state_case(c)
        self.assertFalse(any(x['match_basis']=='book_integer_clue_only' for x in r['observations'][0]['source_leads']))
