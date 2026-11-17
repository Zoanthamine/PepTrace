import unittest
from fusion_bridge import diagnose_fused_case,adapt_dsh_case,peer_assets
from natural_case import parse_text

def case(seq,delta,rules):return {'intended_product':{'sequence':seq,'nterm':'free','cterm':'amide'},
    'reaction':{'description':'fusion regression'},'observations':[{'kind':'delta_mono','value':delta,'tolerance':0.00001,'tolerance_unit':'Da'}],
    'rule_ids':rules,'max_events':1}

class FusionTests(unittest.TestCase):
    def test_canonical_cam_and_acm(self):
        lib,rs=peer_assets();self.assertAlmostEqual(lib.deltas['Carbamidomethyl'].mass,57.021464,places=5)
        self.assertAlmostEqual(lib.deltas['Acm'].mass,71.037114,places=5);self.assertFalse(rs.problems)
    def test_phe_donor_parallel_nalpha_and_lys(self):
        c=case('SYGRKA',147.06841391407,['CH10-R003'])
        c['reaction'].update(stage_tags=['coupling','lys_sidechain_acylation'],acyl_donors=[{'sequence':'F','nterm':'free','cterm':'acid'}])
        r=diagnose_fused_case(c);bonds=[b for h in r['observations'][0]['candidates'] for b in h['component']['state']['bonds']]
        self.assertTrue(any(e.get('site')=='N' for b in bonds for e in b['endpoints']))
        self.assertTrue(any(e.get('site')==5 and e.get('atom')=='Lys_Nepsilon' for b in bonds for e in b['endpoints']))
        self.assertTrue(r['fusion_audit']['whole_product_formula_checks_pass'])
    def test_no_dcc_excludes_dcc_route_but_not_all_trisulfides(self):
        c=case('CGYSCRA',31.9720711744,['CH13-R011','CH13-R014'])
        c['intended_product']['bonds']=[{'kind':'disulfide','sites':[1,5]}]
        c['reaction'].update(stage_tags=['cystine_coupling'],stage_tags_absent=['dcc'])
        o=diagnose_fused_case(c)['observations'][0]
        self.assertTrue(o['candidates']);self.assertTrue(o['condition_conflicts'])
        self.assertTrue(any('dcc' in e['conflicting_conditions'] for h in o['condition_conflicts'] for p in h['paths'] for e in p))
        self.assertFalse(any(e['branch_id']=='CH13-R014:000' for h in o['candidates'] for p in h['paths'] for e in p))
    def test_negative_dcc_and_aliases_without_substring_false_positive(self):
        r=parse_text('序列 CGC，N端游离，C端酰胺，ΔM=31.972071 Da，容差0.001 Da，未加DCC，用甲醇及TFA。')
        self.assertEqual(r['status'],'ready');reaction=r['case']['reaction']
        self.assertIn('dcc',reaction['stage_tags_absent']);self.assertIn('methanol',reaction['stage_tags'])
        self.assertIn('tfa',reaction['stage_tags']);self.assertNotIn('dcc',reaction['stage_tags'])
    def test_no_guessed_adduct_charge_or_second_peak(self):
        for text in ['序列 MELK，N端游离，C端羧酸，m/z=535.1，容差5 ppm。',
            '序列 MELK，N端游离，C端羧酸，m/z=535.1，m/z=536.2，电荷1，[M+H]，容差5 ppm。']:
            self.assertEqual(parse_text(text)['status'],'needs_input')
    def test_legacy_adapter_rejects_implicit_defaults_and_tags(self):
        with self.assertRaises(ValueError):adapt_dsh_case({'sequence':'MELK','observed':{'kind':'mz','value':535}})
        r=adapt_dsh_case({'sequence':'MELK','nterm':'free','cterm':'acid','observed':{'kind':'delta','value':15.994915,'tol_da':0.001},
            'conditions':{'denied_reagents':['DCC'],'reported_reagents':['Methanol']}})
        self.assertEqual(r['reaction']['stage_tags_absent'],['dcc']);self.assertEqual(r['reaction']['stage_tags'],['methanol'])
    def test_peer_nominal_not_promoted_to_high_resolution(self):
        c=case('GAGA',41.5,[]);r=diagnose_fused_case(c)
        self.assertFalse(r['observations'][0]['candidates'])
        self.assertFalse(any(x['mass_grade'].startswith('nominal') for x in r['observations'][0]['dsh_source_review']))
    def test_delete_modified_arg_consumes_its_modification(self):
        c=case('GRA',-198.11167545955,['CH05-R004'])
        c['intended_product']['modifications']=[{'site':2,'name':'Arg乙酰','atom':'guanidino_Nomega','delta_formula':{'C':2,'H':2,'O':1}}]
        r=diagnose_fused_case(c)
        self.assertTrue(r['observations'][0]['candidates'])
    def test_source_or_branches_do_not_overexclude_a_rule(self):
        from knowledge_source import load_runtime
        from source_leads import lookup_delta
        leads=lookup_delta(load_runtime()[0],31.9720711744,1e-5,list('CGYSCRA'),['cystine_coupling','bop'],['dcc'])['source_leads']
        row=next(x for x in leads if x['id']=='CH13-R014')
        self.assertEqual(row['context_status'],'context_supported')
