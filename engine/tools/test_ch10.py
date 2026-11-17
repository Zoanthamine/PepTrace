import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules


class ChapterTenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    @staticmethod
    def target(sequence, mods=(), cterm="acid"):
        return dict(sequence=sequence, nterm="free", cterm=cterm, modifications=list(mods))

    def result(self, target, value, tags=(), donors=None):
        reaction = dict(description="CH10 independent test", stage_tags=list(tags))
        if donors is not None:
            reaction["acyl_donors"] = donors
        return diagnose({"intended_product": target, "reaction": reaction,
                         "observations": [{"kind": "neutral_mono", "value": value,
                                           "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)

    @staticmethod
    def events(result, number):
        return [e for h in result["observations"][0]["matches"] for e in h["events"]
                if e["rule_id"] == f"CH10-R{number:03}"]

    def test_all_repeat_source_rules_preserve_distinct_origin_sites(self):
        parent = mass.calculate_mass(sequence="GAS")
        for number in [1, 2, 6, 7, 8, 9, 11, 12, 15]:
            with self.subTest(rule=number):
                event, = self.events(self.result(self.target("GAS"), parent+mass.calculate_mass(formula="C3H5NO")), number)
                self.assertEqual(event["site"], 2)
                self.assertEqual(event["product_structure"]["sequence"], list("GAAS"))
        event, = self.events(self.result(self.target("GAS"), parent+2*mass.calculate_mass(formula="C3H5NO")), 1)
        self.assertEqual(event["product_structure"]["sequence"], list("GAAAS"))
        # A repeated protected site's downstream phosphate must shift twice and survive.
        target = self.target("GAS", [{"site": 3, "group": "Phospho"}])
        observed = parent+2*mass.calculate_mass(formula="C3H5NO")+mass.calculate_mass(formula="HPO3")
        event, = self.events(self.result(target, observed), 1)
        self.assertEqual(event["product_structure"]["modifications"], [{"site": 5, "group": "Phospho"}])

    def test_fmoc_loss_sources_require_actual_fmoc_and_lys_when_applicable(self):
        observed = mass.calculate_mass(sequence="GKA")
        target = self.target("GKA", [{"site": "N", "group": "Fmoc"}])
        for number in [2, 3, 6, 7, 8, 9, 10, 11]:
            with self.subTest(rule=number):
                event, = self.events(self.result(target, observed), number)
                self.assertEqual(event["delta_formula"], {"C": -15, "H": -10, "O": -2})
                self.assertEqual(event["site"], "N")
                self.assertFalse(self.events(self.result(self.target("GKA"), observed), number))
        self.assertFalse(self.events(self.result(self.target("GAA", [{"site": "N", "group": "Fmoc"}]),
                                                 mass.calculate_mass(sequence="GAA")), 3))

    def test_fragment_self_condensation_excludes_acceptor_and_does_not_add_to_full_target(self):
        # Intended Fmoc-GP-Hyp; product is Fmoc-GP-GP-OH.
        target = self.target(["G", "P", "Hyp"], [{"site": "N", "group": "Fmoc"}])
        donor = self.target("GP", [{"site": "N", "group": "Fmoc"}])
        observed = mass.calculate_mass(formula="C29H32N4O7")
        result = self.result(target, observed, ["pro_hyp_induced_fmoc_deprotection", "segment_condensation"], [donor])
        event, = self.events(result, 4)
        self.assertEqual(event["product_structure"]["sequence"], list("GPGP"))
        self.assertEqual(event["product_structure"]["cterm"], "acid")
        self.assertEqual(event["delta_formula"], {"C": 2, "H": 3, "N": 1})
        self.assertEqual(event["condition"], "compatible")
        missing = self.result(target, observed)
        self.assertIn("CH10-R004", missing["scope"]["component_dependent_rules_skipped"])
        self.assertFalse(self.events(missing, 4))
        self.assertFalse(self.events(self.result(self.target(["G", "P", "A"], [{"site": "N", "group": "Fmoc"}]),
                                                 observed, donors=[donor]), 4))

    def test_pro_and_ctc_gly_sites_are_sequence_limited(self):
        observed = mass.calculate_mass(sequence="PPPA")+mass.calculate_mass(formula="C5H7NO")
        found = self.events(self.result(self.target("PPPA"), observed), 5)
        self.assertEqual({e["site"] for e in found}, {1, 2})
        self.assertTrue(all(e["product_structure"]["sequence"] == list("PPPPA") for e in found))
        parent = mass.calculate_mass(sequence="GAG")
        event, = self.events(self.result(self.target("GAG"), parent+mass.calculate_mass(formula="C2H3NO")), 14)
        self.assertEqual(event["site"], 3)
        self.assertEqual(event["product_structure"]["sequence"], list("GAGG"))
        self.assertFalse(self.events(self.result(self.target("GAG", cterm="amide"),
                                                 parent-0.984015583+mass.calculate_mass(formula="C2H3NO")), 14))

    def test_his_gly_acyl_transfer_and_aminooxy_repeat_have_defined_composition(self):
        parent = mass.calculate_mass(sequence="GHA")
        found = self.events(self.result(self.target("GHA"), parent+mass.calculate_mass(formula="C2H3NO")), 13)
        self.assertEqual(len(found), 2)
        self.assertTrue(all(e["delta_formula"] == {"C": 2, "H": 3, "N": 1, "O": 1} for e in found))
        target = self.target(["Aoa", "G"])
        # 2 x Aoa residue C2H3NO2 + Gly residue C2H3NO + H2O.
        event, = self.events(self.result(target, mass.calculate_mass(formula="C6H11N3O6")), 16)
        self.assertEqual(event["product_structure"]["sequence"], ["Aoa", "Aoa", "G"])
        self.assertEqual(event["delta_formula"], {"C": 2, "H": 3, "N": 1, "O": 2})


if __name__ == "__main__":
    unittest.main()
