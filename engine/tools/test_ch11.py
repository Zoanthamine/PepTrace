import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules, structure


class ChapterElevenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def result(self, target, value, tags=(), features=()):
        return diagnose({"intended_product": target,
                         "reaction": {"description": "CH11 alpha configuration", "stage_tags": list(tags)},
                         "context": {"features_present": list(features)},
                         "observations": [{"kind": "neutral_mono", "value": value,
                                           "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)

    @staticmethod
    def target(seq, **extra):
        return dict(sequence=seq, nterm="free", cterm="acid", modifications=[], **extra)

    @staticmethod
    def events(result, n):
        return [e for h in result["observations"][0]["matches"] for e in h["events"]
                if e["rule_id"] == f"CH11-R{n:03}"]

    def test_all_epimerization_families_preserve_mass_and_report_relative_configuration(self):
        for n, seq, count, mods in [(1, "AGA", 2, []), (2, "AGA", 2, []), (3, "AGA", 2, []),
                                    (5, "AGA", 2, []), (6, "GHA", 2, []), (7, "GCA", 1, []),
                                    (8, "GSA", 1, [{"site": 2, "group": "HexNAc"}]), (9, "GAS", 2, []),
                                    (11, "AGA", 2, []), (12, "GCH", 2, []), (13, "GAS", 2, [])]:
            with self.subTest(rule=n):
                target = self.target(seq)
                target["modifications"] = mods
                observed = mass.calculate_mass(sequence=seq)+(mass.calculate_mass(formula="C8H13NO5") if mods else 0)
                events = self.events(self.result(target, observed), n)
                self.assertEqual(len(events), count)
                self.assertTrue(all(e["delta_formula"] == {} for e in events))
                self.assertTrue(all(e["stereochemistry_edit"]["absolute_assignment"] == "not_assigned" for e in events))
                self.assertTrue(all(seq[e["site"]-1] != "G" for e in events))

    def test_known_ld_inversion_preserves_other_centers_and_round_trips(self):
        target = self.target("ACS", alpha_stereochemistry={"1": "L", "2": "D", "3": "L"})
        event, = self.events(self.result(target, mass.calculate_mass(sequence="ACS")), 7)
        self.assertEqual(event["product_structure"]["alpha_stereochemistry"], {"1": "L", "2": "L", "3": "L"})
        self.assertEqual(event["stereochemistry_edit"]["input_LD"], "D")
        self.assertEqual(event["stereochemistry_edit"]["product_LD"], "L")
        self.assertEqual(structure(event["product_structure"])[1], structure(target)[1])
        self.assertTrue(self.result(target, mass.calculate_mass(sequence="ACS"))["scope"]["sequence_rules_skipped"])
        with self.assertRaisesRegex(ValueError, "chiral alpha"):
            structure(self.target("G", alpha_stereochemistry={"1": "L"}))

    def test_aspartimide_hydrolysis_keeps_asp_and_asn_references_distinct(self):
        for seq, increment, expected in [("GDA", 0, 3), ("GNA", 0.984015583, 4)]:
            with self.subTest(sequence=seq):
                result = self.result(self.target(seq), mass.calculate_mass(sequence=seq)+increment)
                found = self.events(result, 4)
                self.assertEqual(len(found), expected)
                self.assertEqual(sum(e["product_structure"].get("linkage_edit", {}).get("kind") == "isoAsp" for e in found), 2)
                self.assertEqual({e["stereochemistry_edit"]["operation"] for e in found},
                                 {"invert_relative_to_input", "retain_relative_to_input"})
        self.assertFalse(self.events(self.result(self.target("GNP"), mass.calculate_mass(sequence="GNP")+0.984015583), 4))

    def test_c_terminal_cys_ester_piperidine_risk_excludes_amide_and_internal_cys_phase(self):
        for seq, cterm, stereo_feature in [("AGC", "acid", "c_terminal_cys_ester_used:3"),
                                          ("AGC", "amide", "c_terminal_cys_ester_used:3"),
                                          ("ACG", "acid", "c_terminal_cys_ester_used:2")]:
            target = self.target(seq)
            target["cterm"] = cterm
            observed = mass.calculate_mass(sequence=seq)-(0.984015583 if cterm == "amide" else 0)
            events = self.events(self.result(target, observed, ["piperidine_deprotection"], [stereo_feature]), 7)
            supported = [e for e in events if e["condition"] == "compatible"]
            self.assertEqual(len(supported), 1 if seq == "AGC" and cterm == "acid" else 0)

    def test_terminal_aib_racemizes_penultimate_center_and_phg_is_chiral(self):
        target = self.target(["G", "A", "Aib"])
        observed = mass.calculate_mass(formula="C9H17N3O4")
        events = self.events(self.result(target, observed, ["terminal_disubstituted_amino_acid_activation"]), 9)
        supported = [e for e in events if e["condition"] == "compatible"]
        self.assertEqual([e["site"] for e in supported], [2])
        self.assertFalse(any(e["site"] == 3 for e in events))
        phg = self.target(["G", "Phg", "G"])
        event, = self.events(self.result(phg, mass.calculate_mass(formula="C12H15N3O4")), 10)
        self.assertEqual(event["site"], 2)

    def test_control_facts_are_excluded_from_mass_candidates_and_zero_mass_leads(self):
        result = diagnose({"reaction": {"description": "zero mass source lookup"},
                           "observations": [{"kind": "nominal_delta", "value": 0,
                                             "reference": "intended_product"}]}, self.library)
        ids = {r["rule_id"] for r in result["observations"][0]["reference_only_delta_leads"]}
        controls = {f"CH11-R{n:03}" for n in range(14, 35)}
        self.assertFalse(ids & controls)


if __name__ == "__main__":
    unittest.main()
