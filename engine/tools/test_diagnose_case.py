import unittest

from diagnose_case import coordinate, diagnose, load_rules, structure
from mass_ledger import formula_mass


class DiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def run_case(self, sequence, mass, *, cterm="acid", nterm="free",
                 start=None, tags=(), depth=1, observation=None):
        target = {"sequence": sequence, "nterm": nterm, "cterm": cterm}
        case = {"starting_material": start or target, "intended_product": target,
                "reaction": {"description": "synthetic test", "stage_tags": list(tags)},
                "observations": [observation or {"kind": "neutral_mono", "value": mass,
                                                "tolerance": 0.000002, "tolerance_unit": "Da"}],
                "max_events": depth}
        return diagnose(case, self.library)

    def matches(self, *args, **kwargs):
        return self.run_case(*args, **kwargs)["observations"][0]["matches"]

    @staticmethod
    def ids(hits):
        return {e["rule_id"] for h in hits for e in h["events"]}

    def test_target_reference_and_unreacted(self):
        # Free glycine 75.032028, acetylglycine 117.042593.
        self.assertEqual(["starting_material"], [h["kind"] for h in self.matches(
            "G", 75.032028, nterm="acetyl",
            start={"sequence": "G", "nterm": "free", "cterm": "acid"})])
        self.assertIn("intended_product", [h["kind"] for h in self.matches("G", 117.042593, nterm="acetyl")])

    def test_no_asn_or_glutamine(self):
        self.assertEqual([], self.matches("G", 76.016044))

    def test_terminal_and_side_chain_isobaric(self):
        # Free Asn and asparagine amide undergoing side-chain deamidation have equal mass.
        hits = self.matches("N", 132.053492, cterm="amide")
        self.assertEqual({"SUP-001", "SUP-045"}, self.ids(hits))

    def test_nterm_blocked(self):
        # Ac-Gln, less NH3; free N-terminal Gln cyclization does not apply.
        self.assertNotIn("SUP-005", self.ids(self.matches("Q", 171.053158, nterm="acetyl")))

    def test_reagent_gate(self):
        self.assertIn("CH02-R002", self.ids(self.matches("GC", 229.142641, tags=["piperidine"])))
        # A positive tag for another reagent is not evidence that piperidine was never used.
        self.assertIn("CH02-R002", self.ids(self.matches("GC", 229.142641, tags=["4-methylpiperidine"])))

    def test_independent_double_vs_single_site(self):
        hits = self.matches("MM", 312.081364, depth=2)
        self.assertTrue(any([e["rule_id"] for e in h["events"]] == ["CH09-R029"] * 2 for h in hits))
        self.assertFalse(any(len(h["events"]) == 2 for h in self.matches("M", 181.040879, depth=2)))

    def test_ion_coordinates(self):
        for z, ions, expected in [(1, {"H": 1}, 76.039305), (1, {"Na": 1}, 98.021249),
                                  (-1, {"-H": 1}, 74.024752), (2, {"H": 1, "K": 1}, 57.501231)]:
            with self.subTest(ions=ions):
                self.assertTrue(self.matches("G", None, observation={
                    "kind": "mz", "value": expected, "charge": z, "ions": ions,
                    "tolerance": 0.000002, "tolerance_unit": "Da"}))

    def test_missing_or_inconsistent_charge(self):
        observation = {"kind": "mz", "value": 75, "tolerance": 1, "tolerance_unit": "Da", "ions": {"H": 1}}
        with self.assertRaises(KeyError):
            coordinate(observation)
        with self.assertRaises(ValueError):
            coordinate(dict(observation, charge=2))

    def test_known_unreviewed_rules_are_only_leads(self):
        output = self.run_case("G", 99.032028)["observations"][0]
        self.assertFalse(output["matches"])
        self.assertNotIn("CH07-R031", {item["rule_id"] for item in output["reference_only_delta_leads"]})

    def test_named_groups_match_supplier_building_block_formulas(self):
        for residue, group, expected in [
            ("R", "Pmc", {"C": 35, "H": 42, "N": 4, "O": 7, "S": 1}),
            ("K", "Dde", {"C": 31, "H": 36, "N": 2, "O": 6}),
        ]:
            with self.subTest(group=group):
                _, formula, _ = structure({"sequence": residue, "nterm": "free", "cterm": "acid",
                                            "modifications": [{"site": "N", "group": "Fmoc"},
                                                              {"site": 1, "group": group}]})
                self.assertEqual(expected, formula)
        self.assertAlmostEqual(699.293799, formula_mass(structure({
            "sequence": "GCK", "nterm": "free", "cterm": "acid",
            "modifications": [{"site": "N", "group": "Fmoc"}, {"site": 2, "group": "Acm"},
                              {"site": 3, "group": "Boc"}]})[1]), places=6)

    def test_partial_deprotection_requires_explicit_starting_group(self):
        # Gly-Lys free acid + Fmoc only, with Boc already removed.
        start = {"sequence": "GK", "nterm": "free", "cterm": "acid",
                 "modifications": [{"site": "N", "group": "Fmoc"}, {"site": 2, "group": "Boc"}]}
        hits = self.matches("GK", 425.195071, start=start, tags=["deprotection"])
        self.assertEqual({"SUP-014", "CH14-R002"}, self.ids(hits))
        self.assertEqual(set(), self.ids(self.matches("GK", 425.195071, tags=["deprotection"])))


if __name__ == "__main__":
    unittest.main()
