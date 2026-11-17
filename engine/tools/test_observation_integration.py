import unittest
from diagnose_case import diagnose

class ObservationIntegrationTests(unittest.TestCase):
    def case(self,observation):
        return diagnose({'intended_product':{'sequence':'G','nterm':'free','cterm':'acid'},
                         'reaction':{'description':'coordinate integration'},'observations':[observation]})

    def test_average_neutral_mass_uses_average_composition(self):
        # C2H5NO2: 2*12.011 + 5*1.008 +14.007+2*15.999=75.067.
        r=self.case({'kind':'neutral_average','value':75.067,'tolerance':0.000001,'tolerance_unit':'Da'})
        h=next(h for h in r['observations'][0]['matches'] if h['kind']=='intended_product')
        self.assertAlmostEqual(h['predicted'],75.067,places=7)
        self.assertAlmostEqual(h['neutral_mono'],75.03202840472,places=7)
        self.assertEqual(h['observation_mass_basis'],'average')

    def test_declared_c13_shift_matches_the_same_molecule(self):
        r=self.case({'kind':'mz','value':77.042659706341,'charge':1,'ions':{'H':1},
                     'isotope_shift':1.003354835,'tolerance':0.000001,'tolerance_unit':'Da'})
        self.assertTrue(any(h['kind']=='intended_product' for h in r['observations'][0]['matches']))

    def test_mono_mass_does_not_silently_accept_average_label(self):
        with self.assertRaisesRegex(ValueError,'disagree'):
            self.case({'kind':'neutral_mono','mass_basis':'average','value':75.067,'tolerance':0.1,'tolerance_unit':'Da'})
