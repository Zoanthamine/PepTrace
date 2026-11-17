import json,math,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from _experimental import spps_review,hplc_review,cyclic_review
from _ms1 import extract_scans,mzml_scans
from _protocols import normalize_snapshot,fetch_snapshot,validate_fetch_url,_NoRedirect

class SPPSTests(unittest.TestCase):
    def test_sequence_and_unknown_process_are_separate(self):
        r=spps_review({'peptide':{'sequence':'QDGCM'}})
        self.assertIsNone(r['predicted_success_probability']);self.assertIsNone(r['reference_element_mass'])
        self.assertEqual(r['length'],5);self.assertTrue(r['missing_process_fields']);self.assertEqual(r['missing_terminal_fields'],['nterm','cterm'])
        self.assertTrue(any(x['name']=='Asp-X基序' and x['sites']==[2,3] for x in r['findings']))
    def test_mass_uses_the_canonical_terminal_state(self):
        from pyteomics import mass
        r=spps_review({'peptide':{'sequence':'GC','nterm':'free','cterm':'amide'},'process':{'loading_mmol_g':.3}})
        self.assertAlmostEqual(r['reference_element_mass'],mass.calculate_mass(sequence='GC')-.98401558291,places=6)
    def test_bad_sequence_or_loading_is_rejected(self):
        for raw in [{'peptide':{'sequence':'GX'}},{'peptide':{'sequence':'GC'},'process':{'loading_mmol_g':-1}}]:
            with self.assertRaises(ValueError):spps_review(raw)
    def test_asp_pro_is_not_given_the_generic_next_nh_aspartimide_flag(self):
        r=spps_review({'peptide':{'sequence':'ADPG'}})
        self.assertTrue(any(x['name']=='Asp-Pro酸解位点' for x in r['findings']))
        self.assertFalse(any(x['name']=='Asp-X基序' for x in r['findings']))

def method():return {'purpose':'LC-MS','column':{'length_mm':150,'id_mm':4.6,'pH_range':[2,8]},
    'flow_ml_min':1,'gradient':[{'time_min':0,'percent_B':5},{'time_min':30,'percent_B':65}],
    'mobile_phase':{'organic':'methanol','additives':['TFA'],'pH':2.5},'transfer_to':{'length_mm':100,'id_mm':2.1}}
class HPLCTests(unittest.TestCase):
    def test_gradient_scaling_preserves_column_volume_basis(self):
        r=hplc_review(method());transfer=r['transfer']
        self.assertAlmostEqual(r['geometric_column_volume_ml'],math.pi*2.3**2*.15)
        self.assertAlmostEqual(transfer['suggested_flow_ml_min'],(2.1/4.6)**2)
        self.assertAlmostEqual(transfer['gradient'][-1]['time_min'],20);self.assertEqual(r['gradient_slopes_percent_B_per_min'],[2.0])
    def test_conditions_are_kept_as_findings_not_retention_predictions(self):
        r=hplc_review(method());names={x['name'] for x in r['findings']}
        self.assertTrue({'TFA/ESI兼容性','酸/甲醇甲酯化线索'}<=names)
    def test_invalid_gradient_and_column_limits(self):
        raw=method();raw['gradient'][1]['time_min']=0
        with self.assertRaises(ValueError):hplc_review(raw)
        raw=method();raw['mobile_phase']['pH']=9
        self.assertTrue(any(x['name']=='超出声明的柱pH范围' for x in hplc_review(raw)['findings']))

RING='[CH2:1]1[NH:2][C:3](=[O:4])[CH2:5][NH:6][C:7]1=[O:8]'
class CyclicTests(unittest.TestCase):
    def test_formal_charge_does_not_claim_permanent_intrinsic_charge(self):
        r=cyclic_review({'smiles':'C[NH3+]','expected_formal_charge':1})
        self.assertEqual(r['formal_charge'],1);self.assertTrue(r['expected_charge_matches']);self.assertNotIn('intrinsic_charge',r)
    def test_actual_amide_ring_and_acid_precursor_water_balance(self):
        r=cyclic_review({'smiles':RING,'closure':{'kind':'amide','atom_maps':[2,3]},'expected_formula':'C4H6N2O2',
            'peptide_reference':{'sequence':'GG','nterm':'free','cterm':'acid'},'reference_relation':'head_to_tail'})
        self.assertTrue(r['closure_check']['pass']);self.assertTrue(r['expected_formula_matches']);self.assertTrue(r['reference_relation']['matches_actual_formula'])
    def test_disulfide_is_not_an_amide_closure(self):
        smiles='[CH2:1]1[CH2:2][S:3][S:4][CH2:5]1'
        self.assertTrue(cyclic_review({'smiles':smiles,'closure':{'kind':'disulfide','atom_maps':[3,4]}})['closure_check']['pass'])
        self.assertFalse(cyclic_review({'smiles':smiles,'closure':{'kind':'amide','atom_maps':[3,4]}})['closure_check']['pass'])
    def test_wrong_formula_or_amidated_reference_is_not_accepted(self):
        self.assertFalse(cyclic_review({'smiles':RING,'expected_formula':'C4H8N2O3'})['expected_formula_matches'])
        with self.assertRaises(ValueError):cyclic_review({'smiles':RING,'reference_relation':'head_to_tail','peptide_reference':{'sequence':'GG','nterm':'free','cterm':'amide'}})

def targets():return {'targets':[{'id':'x','mz':500.0,'tolerance':10,'tolerance_unit':'ppm'}]}
def scans():return [{'ms_level':1,'centroided':True,'rt_seconds':t,'peaks':[[500.001,i],[600.,1.]]} for t,i in [(0.,0.),(10.,10.),(20.,0.)]]
class MS1Tests(unittest.TestCase):
    def test_xic_tic_and_area_against_independent_triangle(self):
        r=extract_scans(scans(),targets());self.assertEqual([x['TIC'] for x in r['ms1_rows']],[1.,11.,1.]);self.assertEqual(r['area_intensity_seconds']['x'],100.)
    def test_overlapping_targets_are_not_exclusive_assignments(self):
        spec=targets();spec['targets'].append({'id':'y','mz':500.002,'tolerance':.01,'tolerance_unit':'Da'})
        r=extract_scans(scans(),spec);self.assertEqual(r['overlapping_peak_target_assignments'],3)
        self.assertEqual(r['area_intensity_seconds']['x'],r['area_intensity_seconds']['y'])
    def test_profile_negative_intensity_and_decreasing_rt_rejected(self):
        cases=[scans(),scans(),list(reversed(scans()))];cases[0][0]['centroided']=False;cases[1][0]['peaks'][0][1]=-1
        for s in cases:
            with self.assertRaises(ValueError):extract_scans(s,targets())
    def test_indexed_mzml_reads_match_saved_centroid_scans(self):
        import pyopenms as ms
        work=ROOT/'validation/work';work.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=work) as folder:
            path=Path(folder)/'ms1.mzML';exp=ms.MSExperiment()
            for i,row in enumerate(scans()):
                s=ms.MSSpectrum();s.setMSLevel(1);s.setRT(row['rt_seconds']);s.setNativeID('scan='+str(i+1));s.setType(ms.SpectrumSettings.SpectrumType.CENTROID)
                s.set_peaks(([p[0] for p in row['peaks']],[p[1] for p in row['peaks']]));exp.addSpectrum(s)
            ms.MzMLFile().store(str(path),exp)
            self.assertEqual(extract_scans(mzml_scans(path),targets())['area_intensity_seconds']['x'],100.)

class ProtocolTests(unittest.TestCase):
    def test_version_guid_order_and_html_remain_explicit(self):
        raw={'protocol':{'title':'Example','doi':'10.17504/protocols.io.example/v2','version_uri':'example/v2',
            'authors':[{'name':'A'}],'steps':[{'guid':'a','previous_guid':None,'step':'<p>Wash resin</p><script>execute()</script>'},
                {'guid':'b','previous_guid':'a','step':'Measure'}]}}
        r=normalize_snapshot(raw);self.assertEqual(r['steps'][0]['text'],'Wash resin');self.assertTrue(r['guid_links_consistent_with_source_order'])
        self.assertEqual(r['version_uri'],'example/v2');self.assertEqual(r['authors'],[{'name':'A'}]);self.assertTrue(r['untrusted_content'])
    def test_payload_and_unsupported_content_are_reported(self):
        r=normalize_snapshot({'payload':{'title':'x','steps':[{'content':{'unknown':'shape'}}]}})
        self.assertEqual(r['unsupported_step_content_count'],1);self.assertIn('version_uri',r['missing_provenance'])
    def test_default_fetch_and_missing_token_do_not_make_network_calls(self):
        url='https://www.protocols.io/api/v4/protocols/123'
        self.assertFalse(fetch_snapshot(url)['network_performed'])
        with patch.dict('os.environ',{},clear=True):self.assertEqual(fetch_snapshot(url,True)['status'],'needs_access_token')
    def test_unknown_host_and_redirect_are_rejected(self):
        for url in ['https://example.org/api/v4/protocols/123','http://www.protocols.io/api/v4/protocols/123','https://user:secret@www.protocols.io/api/v4/protocols/123']:
            with self.assertRaises(ValueError):validate_fetch_url(url)
        with self.assertRaises(ValueError):_NoRedirect().redirect_request(None,None,302,'',{},'https://example.org')
