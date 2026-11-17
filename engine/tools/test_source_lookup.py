import unittest
from knowledge_source import load_runtime
from source_leads import lookup_delta
class SourceLookupTests(unittest.TestCase):
    def test_missing_termini_and_precision_can_return_integer_source_clues(self):
        r=lookup_delta(load_runtime()[0],15.99,None,list('MELK'))
        self.assertFalse(r['state_search_performed'])
        self.assertTrue(any(x['id']=='CH09-R029' for x in r['source_leads']))
        self.assertTrue(all(x['classification']=='source_lead_only' for x in r['source_leads']))
        self.assertIn('未判断高分辨',r['precision_note'])
