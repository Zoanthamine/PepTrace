import unittest
from natural_case import parse_text
from state_diagnosis import diagnose_state_case
from state_report import render_state_report

class NaturalStateCaseTests(unittest.TestCase):
    def test_explicit_text_to_traceable_report(self):
        r=parse_text('序列 MELK，N端游离，C端羧酸，实测比理论多15.994915 Da，容差0.001 Da，用DMSO处理过。')
        self.assertEqual(r['status'],'ready')
        c=r['case'];c['rule_ids']=['CH09-R029'];c['max_events']=1
        d=diagnose_state_case(c);md=render_state_report(d)
        self.assertTrue(d['observations'][0]['candidates'])
        self.assertIn('CH09-R029',md);self.assertIn('223',md)
        self.assertEqual(d['case']['input_text'],r['case']['input_text'])
    def test_missing_inputs_are_not_guessed(self):
        r=parse_text('序列 MELK，比理论多15.99 Da，用DMSO处理。')
        self.assertEqual(r['status'],'needs_input')
        self.assertIn('N端状态',r['missing']);self.assertIn('实际 Da／ppm 容差',r['missing'])
    def test_denied_condition_and_negative_delta(self):
        r=parse_text('序列 GC，N端游离，C端酸，质量差-33.987721 Da，容差0.001 Da，未使用哌啶。')
        self.assertEqual(r['status'],'ready');self.assertEqual(r['case']['observations'][0]['value'],-33.987721)
        self.assertIn('piperidine',r['case']['reaction']['stage_tags_absent'])
    def test_mz_charge_and_adduct_are_retained(self):
        r=parse_text('序列 G，N端游离，C端酸，高分辨 m/z=76.039305，电荷1，质子化[M+H]，容差5 ppm。')
        self.assertEqual(r['status'],'ready');self.assertEqual(r['case']['observations'][0]['ions'],{'H':1})
