import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules


class ChapterEightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def result(self, target, value, tags=(), donors=None):
        reaction = {"description": "CH08 independent formula tests", "stage_tags": list(tags)}
        if donors is not None:
            reaction["acyl_donors"] = donors
        return diagnose({"intended_product": target, "reaction": reaction,
                         "observations": [{"kind": "neutral_mono", "value": value,
                                           "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)

    @staticmethod
    def target(seq, mods=(), cterm="acid", nterm="free"):
        return {"sequence": seq, "nterm": nterm, "cterm": cterm, "modifications": list(mods)}

    @staticmethod
    def events(result, number):
        return [e for h in result["observations"][0]["matches"] for e in h["events"]
                if e["rule_id"] == f"CH08-R{number:03}"]

    def test_local_rules_independent_deltas_and_residue_selectivity(self):
        parent = mass.calculate_mass(sequence="GDEST")
        for num, delta, sites in [
            (1, "CH2", {2, 3}), (4, "CH2", {2, 3}), (6, "CH2", {2, 3}),
            (8, "C3H4", {4, 5}), (12, "CO", {4, 5}), (13, "C2H2O", {4, 5}),
            (14, "C2H-1F3O", {4, 5}), (16, "H-2O-1", {4}), (17, "H-2O-1", {5}),
            (18, "C5H9NO-1", {4}), (19, "C5H9NO-1", {5}),
            (26, "C-1H-2O-1", {4}), (27, "C-2H-4O-1", {5}),
        ]:
            with self.subTest(rule=num):
                found = self.events(self.result(self.target("GDEST"), parent + mass.calculate_mass(formula=delta)), num)
                self.assertEqual({e["site"] for e in found}, sites)
                self.assertFalse(self.events(self.result(self.target("GKA"), mass.calculate_mass(sequence="GKA") + mass.calculate_mass(formula=delta)), num))

    def test_protected_ester_references_are_distinct(self):
        parent = mass.calculate_mass(sequence="GDE")
        target = self.target("GDE", [{"site": 2, "group": "Bzl"}])
        methyl, = self.events(self.result(target, parent + mass.calculate_mass(formula="CH2"),
                                          ["methanol", "tmah", "chloromethyl_resin_loading"]), 1)
        self.assertEqual(methyl["delta_formula"], {"C": -6, "H": -4})
        tert, = self.events(self.result(target, parent + mass.calculate_mass(formula="C4H8"),
                                        ["tmah", "tert_butanol"]), 2)
        self.assertEqual(tert["delta_formula"], {"C": -3, "H": 2})
        self.assertEqual(tert["condition"], "compatible")
        self.assertFalse(self.events(self.result(self.target("GDE"), parent + mass.calculate_mass(formula="C4H8")), 2))

    def test_carboxyl_end_cannot_be_esterified_when_amidated(self):
        parent = mass.calculate_mass(sequence="GDE")
        for n, formula in [(3, "C4H8"), (5, "CH2"), (7, "CH2")]:
            found, = self.events(self.result(self.target("GDE"), parent + mass.calculate_mass(formula=formula)), n)
            self.assertEqual(found["site"], "C")
            self.assertFalse(self.events(self.result(self.target("GDE", cterm="amide"),
                                                      parent - 0.984015583 + mass.calculate_mass(formula=formula)), n))

    def test_o_alloc_and_o_allyl_are_not_same_group(self):
        parent = mass.calculate_mass(sequence="GST")
        target = self.target("GST", [{"site": 2, "group": "Alloc"}])
        event, = self.events(self.result(target, parent + mass.calculate_mass(formula="C3H4"), ["alloc_deprotection", "pd0"]), 8)
        self.assertEqual(event["delta_formula"], {"C": -1, "O": -2})
        self.assertEqual(event["site"], 2)

    def test_explicit_acyl_donor_preserves_its_protection_and_missing_donor_is_visible(self):
        target = self.target("GSG")
        donor = self.target("V", [{"site": "N", "group": "Fmoc"}])
        observed = mass.calculate_mass(formula="C27H32N4O8")
        tags = ["coupling", "minimal_protection"]
        result = self.result(target, observed, tags, [donor])
        event, = self.events(result, 9)
        self.assertEqual(event["product_structure"]["chains"][1]["modifications"], donor["modifications"])
        self.assertEqual(event["condition"], "compatible")
        self.assertEqual(event["site"], [[1, 2, "sidechain_O"], [2, "C", "carbonyl_C"]])
        missing = self.result(target, observed, tags)
        self.assertIn("CH08-R009", missing["scope"]["component_dependent_rules_skipped"])
        self.assertFalse(self.events(missing, 9))
        with self.assertRaisesRegex(ValueError, "free carboxylic acid"):
            self.result(target, observed, tags, [self.target("V", cterm="amide")])

    def test_n_to_o_transfer_of_three_defined_acyl_groups(self):
        parent = mass.calculate_mass(sequence="SG")
        for group, formula in [("Ac", "C2H2O"), ("Formyl", "CO"), ("Palmitoyl", "C16H30O")]:
            with self.subTest(group=group):
                target = self.target("SG", [{"site": "N", "group": group}])
                event, = self.events(self.result(target, parent + mass.calculate_mass(formula=formula), ["acid_exposure"]), 15)
                self.assertEqual(event["delta_formula"], {})
                self.assertEqual(event["product_structure"]["nterm"], "free")
                self.assertEqual(event["product_structure"]["modifications"][0]["site"], 1)

    def test_oxazolidone_keeps_urethane_carbonyl_and_changes_reference_mass(self):
        observed = mass.calculate_mass(sequence="SG") + mass.calculate_mass(formula="CH-2O")
        for group, delta in [("Moz", {"C": -8, "H": -10, "O": -2}),
                             ("Cbz", {"C": -7, "H": -8, "O": -1}),
                             ("Boc", {"C": -4, "H": -10, "O": -1})]:
            with self.subTest(group=group):
                target = self.target("SG", [{"site": "N", "group": group}])
                event, = self.events(self.result(target, observed, ["base_exposure"]), 21)
                self.assertEqual(event["delta_formula"], delta)
                self.assertEqual(event["product_structure"]["topology"], "local_cycle")
                self.assertEqual(event["site"], ["N", 1])
                blocked = self.target("SG", [{"site": "N", "group": group}, {"site": 1, "group": "tBu"}])
                self.assertFalse(self.events(self.result(blocked, observed + mass.calculate_mass(formula="C4H8")), 21))
        self.assertFalse(self.events(self.result(self.target("SG"), observed), 21))

    def test_backbone_migration_and_two_cyclodehydration_states(self):
        parent = mass.calculate_mass(sequence="GSG")
        for n, delta in [(22, ""), (23, "H-2O-1"), (24, "H-4O-1")]:
            observed = parent + (mass.calculate_mass(formula=delta) if delta else 0)
            event, = self.events(self.result(self.target("GSG"), observed), n)
            self.assertEqual(event["delta_formula"], {} if not delta else {"H": -2 if n == 23 else -4, "O": -1})
            self.assertFalse(self.events(self.result(self.target("SG"), mass.calculate_mass(sequence="SG") + observed - parent), n))

    def test_cys_and_dap_heterocycles_preserve_sulfur_or_nitrogen(self):
        for sequence, parent in [("GCG", mass.calculate_mass(sequence="GCG")),
                                 (["G", "Dap", "G"], mass.calculate_mass(formula="C7H14N4O4"))]:
            for delta, cycle in [("H-2O-1", "thiazoline" if sequence == "GCG" else "imidazoline"),
                                 ("H-4O-1", "thiazole" if sequence == "GCG" else "imidazole")]:
                with self.subTest(sequence=sequence, cycle=cycle):
                    event, = self.events(self.result(self.target(sequence), parent + mass.calculate_mass(formula=delta)), 25)
                    self.assertEqual(event["product_structure"]["cycle"]["kind"], cycle)

    def test_repeat_and_delete_templates_preserve_independent_sites(self):
        parent = mass.calculate_mass(sequence="GSTG")
        event, = self.events(self.result(self.target("GSTG"), parent + mass.calculate_mass(formula="C3H5NO2")), 11)
        self.assertEqual(event["site"], 2)
        self.assertEqual(event["product_structure"]["sequence"], list("GSSTG"))
        event, = self.events(self.result(self.target("GSTG"), parent - mass.calculate_mass(formula="C4H7NO2")), 20)
        self.assertEqual(event["site"], 3)
        self.assertEqual(event["product_structure"]["sequence"], list("GSG"))

    def test_elimination_from_protected_ser_recalculates_actual_leaving_group(self):
        parent = mass.calculate_mass(sequence="GSG")
        observed = parent - mass.calculate_mass(formula="H2O")
        for group, increment in [("Ac", "C2H2O"), ("Formyl", "CO"), ("Boc", "C5H8O2"), ("Phospho", "HPO3")]:
            with self.subTest(group=group):
                target = self.target("GSG", [{"site": 2, "group": group}])
                event, = self.events(self.result(target, observed, ["base_exposure"]), 16)
                self.assertAlmostEqual(mass.calculate_mass(composition=event["delta_formula"]),
                                       -mass.calculate_mass(formula="H2O")-mass.calculate_mass(formula=increment), places=7)
                self.assertEqual(event["product_structure"]["sequence"], ["G", "Dha", "G"])
                self.assertEqual(event["product_structure"]["modifications"], [])


if __name__ == "__main__":
    unittest.main()
