import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules


class ChapterTwelveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def result(self, sequence, mods, value):
        return diagnose({"intended_product": {"sequence": sequence, "nterm": "free", "cterm": "acid", "modifications": mods},
                         "reaction": {"description": "CH12 phosphate chemistry"},
                         "observations": [{"kind": "neutral_mono", "value": value, "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)

    @staticmethod
    def events(result, n):
        return [e for h in result["observations"][0]["matches"] for e in h["events"] if e["rule_id"] == f"CH12-R{n:03}"]

    def test_h_phosphonate_uses_phosphate_reference_and_keeps_named_state(self):
        for residue in "STY":
            sequence = "G"+residue+"G"
            observed = mass.calculate_mass(sequence=sequence)+mass.calculate_mass(formula="HPO2")
            event, = self.events(self.result(sequence, [{"site": 2, "group": "Phospho"}], observed), 1)
            self.assertEqual(event["delta_formula"], {"O": -1})
            self.assertEqual(event["product_structure"]["modifications"][0]["group"], "HPhosphonate")
            self.assertFalse(self.events(self.result(sequence, [], observed), 1))

    def test_pyrophosphate_preserves_two_parent_chains_and_only_phosphorylated_partner_sites(self):
        sequence = "GYYG"
        parent = mass.calculate_mass(sequence=sequence)+2*mass.calculate_mass(formula="HPO3")
        observed = 2*parent-mass.calculate_mass(formula="H2O")
        found = self.events(self.result(sequence, [{"site": 2, "group": "Phospho"}, {"site": 3, "group": "Phospho"}], observed), 2)
        self.assertEqual(len(found), 3)
        self.assertTrue(all(len(e["product_structure"]["chains"]) == 2 for e in found))
        self.assertTrue(all(e["product_structure"]["bonds"][0]["kind"] == "P_O_P_pyrophosphate" for e in found))
        one = mass.calculate_mass(sequence=sequence)+mass.calculate_mass(formula="HPO3")
        found = self.events(self.result(sequence, [{"site": 2, "group": "Phospho"}],
                                        2*one-mass.calculate_mass(formula="H2O")), 2)
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0]["site"], [[1, 2, "phosphate_P"], [2, 2, "phosphate_P"]])
        self.assertFalse(self.events(self.result(sequence, [{"site": 2, "group": "Phospho"}],
                                                 one-mass.calculate_mass(formula="H2O")), 2))

    def test_missing_phosphotyrosine_removes_its_actual_group_and_handles_first_synthesis_position(self):
        for sequence, position, product in [("GYG", 2, "GG"), ("YG", 1, "G")]:
            event, = self.events(self.result(sequence, [{"site": position, "group": "Phospho"}],
                                             mass.calculate_mass(sequence=product)), 3)
            self.assertEqual(event["product_structure"]["sequence"], list(product))
            self.assertEqual(event["product_structure"]["modifications"], [])
            self.assertEqual(event["delta_formula"], {"C": -9, "H": -10, "N": -1, "O": -5, "P": -1})
            self.assertEqual(event["origin"], "inference_from_book")
        self.assertFalse(self.events(self.result("GY", [{"site": 2, "group": "Phospho"}],
                                                 mass.calculate_mass(sequence="G")), 3))


if __name__ == "__main__":
    unittest.main()
