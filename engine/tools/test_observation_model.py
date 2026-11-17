"""Independent physical-coordinate examples; integration with the candidate engine remains R03 work."""
import unittest
from observation_model import formula_mass, hypotheses, parse_coordinate, peak_pair

class ObservationModelTests(unittest.TestCase):
    def test_ppm_boundary_uses_theoretical_coordinate_consistently(self):
        coordinate=parse_coordinate({'kind':'neutral_mono','value':100.001000005,'tolerance':10,'tolerance_unit':'ppm'})
        self.assertFalse(coordinate.match(100)['match'])
    def test_average_is_not_mono_rounding(self):
        self.assertAlmostEqual(formula_mass({'H':2,'O':1},'average'),18.015,places=9)
        self.assertAlmostEqual(formula_mass({'H':2,'O':1}),18.01056468403,places=9)

    def test_fixed_positive_and_negative_charge_electron_balance(self):
        p=parse_coordinate({'kind':'mz','value':999.999451420091,'tolerance':0.000001,'tolerance_unit':'Da','charge':1,'intrinsic_charge':1,'ions':{}})
        self.assertAlmostEqual(p.predict(1000),999.999451420091,places=9)
        n=parse_coordinate({'kind':'mz','value':1000.000548579909,'tolerance':0.000001,'tolerance_unit':'Da','charge':-1,'intrinsic_charge':-1,'ions':{}})
        self.assertAlmostEqual(n.elemental_reference(),1000,places=9)

    def test_double_protonated_and_isotope_roundtrip(self):
        x=parse_coordinate({'kind':'mz','value':501.508953884121,'tolerance':5,'tolerance_unit':'ppm',
                            'charge':2,'ions':{'H':2},'isotope_shift':1.003354835})
        self.assertAlmostEqual(x.elemental_reference(),1000,places=9)
        self.assertAlmostEqual(x.tolerance_da,501.508953884121*5e-6,places=9)

    def test_different_charge_peaks_and_error_propagation(self):
        a=parse_coordinate({'kind':'mz','value':501.007276466621,'tolerance':0.01,'tolerance_unit':'Da','charge':2,'ions':{'H':2}})
        b=parse_coordinate({'kind':'mz','value':1017.002191086191,'tolerance':0.02,'tolerance_unit':'Da','charge':1,'ions':{'H':1}})
        result=peak_pair(a,b)
        self.assertAlmostEqual(result['delta_element_mass'],15.99491461957,places=8)
        self.assertAlmostEqual(result['tolerance_da'],0.04,places=9)
        self.assertAlmostEqual(peak_pair(a,b,'independent_gaussian')['tolerance_da'],0.0282842712474619,places=9)

    def test_declared_ion_hypotheses_are_kept_separately(self):
        options=hypotheses({'kind':'mz','value':500,'tolerance':0.1,'tolerance_unit':'Da',
            'ion_hypotheses':[{'charge':1,'ions':{'H':1}},{'charge':1,'ions':{'Na':1}},{'charge':2,'ions':{'H':2}}]})
        self.assertEqual(len(options),3)
        self.assertEqual(len({x.elemental_reference() for x in options}),3)

    def test_invalid_charge_basis_and_fields_are_rejected(self):
        for raw in [
            {'kind':'mz','value':500,'charge':2,'ions':{'H':1}},
            {'kind':'neutral_mono','value':500,'mass_basis':'average'},
            {'kind':'neutral_average','value':500,'isotope_shift':1.003},
            {'kind':'neutral_mono','value':float('nan')},
            {'kind':'neutral_mono','value':500,'typo':True},
        ]:
            with self.subTest(raw=raw),self.assertRaises(ValueError):
                parse_coordinate({'tolerance':0.01,'tolerance_unit':'Da',**raw})
