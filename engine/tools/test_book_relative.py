"""Regression cases from Yang (2016), independent counterexamples, and input boundaries."""
import unittest

from diagnose_case import diagnose, load_rules


def make_case(value, tags=(), context=None, *, exact=False):
    observation = {"kind": "nominal_delta", "value": value, "reference": "intended_product"}
    if exact:
        observation.update(kind="delta_mono", tolerance_unit="Da", tolerance=0.000002)
    return {"reaction": {"description": "Book regression or counterexample", "stage_tags": list(tags)},
            "context": context or {}, "observations": [observation]}


class RelativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def run_case(self, *a, **kw):
        return diagnose(make_case(*a, **kw), self.library)["observations"][0]

    def test_book_ch14_unknown_structure(self):
        result = self.run_case(27, ["acetonitrile_stress", "alkaline_exposure"], {"features_present": ["ketone"]})
        lead = result["reference_only_delta_leads"][0]
        self.assertEqual("CH14-R012", lead["rule_id"])
        self.assertEqual("context_supported", lead["context_status"])
        self.assertEqual("unconfirmed_hypothesis", lead["structure_status"])
        self.assertIsNone(lead["delta_formula"])
        self.assertIsNone(lead["theoretical_delta_mono"])
        self.assertIsNone(lead["error_da"])
        self.assertFalse(result["matches"])

    def test_approximate_window_preserves_book_only_unknown(self):
        case = make_case(27.2, ["acetonitrile_stress"], {"features_present": ["ketone"]}, exact=True)
        case["observations"][0].update(kind="delta_approx", tolerance=0.3)
        lookup = diagnose(case, self.library)["observations"][0]
        lead = next(row for row in lookup["reference_only_delta_leads"] if row["rule_id"] == "CH14-R012")
        self.assertIsNone(lead["error_da"])
        self.assertEqual(lead["mass_basis"], "book_reported_integer_in_window")

    def test_book_ch03_partial_sequence(self):
        result = self.run_case(172, ["acid_global_deprotection", "edt"], {"residues_present": ["W"]})
        lead = result["reference_only_delta_leads"][0]
        self.assertEqual("CH03-R009", lead["rule_id"])
        self.assertEqual("context_supported", lead["context_status"])
        self.assertAlmostEqual(171.962827, lead["theoretical_delta_mono"], places=6)

    def test_book_ch02_terminal_ambiguity(self):
        result = self.run_case(51, ["base_deprotection", "piperidine"], {"sequence": "DEQEALNSC"})
        lead = result["reference_only_delta_leads"][0]
        self.assertEqual("CH02-R002", lead["rule_id"])
        self.assertAlmostEqual(51.101428, lead["theoretical_delta_mono"], places=6)
        self.assertFalse(result["matches"])

    def test_ketone_unknown(self):
        lead = self.run_case(27, ["acetonitrile_stress"])["reference_only_delta_leads"][0]
        self.assertEqual("context_partial", lead["context_status"])
        self.assertIn("ketone", lead["missing"])

    def test_ketone_explicitly_absent(self):
        lead = self.run_case(27, ["acetonitrile_stress"], {"features_absent": ["ketone"]})["reference_only_delta_leads"][0]
        self.assertEqual("context_conflict", lead["context_status"])

    def test_unconfirmed_nominal_is_not_exact(self):
        self.assertNotIn("CH14-R012", {r["rule_id"] for r in self.run_case(
            27.010899, ["acetonitrile_stress"], {"features_present": ["ketone"]}, exact=True
        )["reference_only_delta_leads"]})

    def test_nominal_27_does_not_round_formyl_down(self):
        self.assertNotIn("SUP-020", {r["rule_id"] for r in self.run_case(27)["reference_only_delta_leads"]})
        self.assertIn("SUP-020", {r["rule_id"] for r in self.run_case(28)["reference_only_delta_leads"]})

    def test_full_sequence_without_cys(self):
        lead = next(r for r in self.run_case(51, ["piperidine"], {"sequence": "MKS"})[
            "reference_only_delta_leads"] if r["rule_id"] == "CH02-R002")
        self.assertEqual("context_conflict", lead["context_status"])

    def test_edt_explicitly_absent(self):
        lead = self.run_case(172, ["acid_global_deprotection"],
                             {"residues_present": ["W"], "stage_tags_absent": ["edt"]})["reference_only_delta_leads"][0]
        self.assertIn("edt", lead["conflicts"])

    def test_negative_book_delta(self):
        lead = self.run_case(-16, ["phosphitylation"],
                             {"features_present": ["phosphopeptide_target"]})["reference_only_delta_leads"][0]
        self.assertEqual("CH12-R001", lead["rule_id"])

    def test_book_proline_insertion_reference(self):
        lead = self.run_case(97, ["fmoc_spps"], {"sequence": "PPPPPPA"})["reference_only_delta_leads"][0]
        self.assertEqual("CH10-R005", lead["rule_id"])

    def test_appendix_trp_separate_from_his_supplement(self):
        result = self.run_case(44)
        appendix = next(r for r in result["appendix_only_leads"] if r["appendix_id"] == "AI-056")
        self.assertIn("Trp", appendix["name"])
        self.assertAlmostEqual(43.989829, appendix["theoretical_delta_mono"], places=6)
        self.assertIn("SUP-053", {r["rule_id"] for r in result["reference_only_delta_leads"]})

    def test_fixed_charge_is_not_neutral(self):
        result = self.run_case(117)
        item = next(r for r in result["appendix_only_leads"] if r["appendix_id"] == "AI-072")
        self.assertIsNone(item["theoretical_delta_mono"])

    def test_unknown_reference_is_rejected(self):
        case = make_case(51)
        case["observations"][0]["reference"] = "starting_material"
        with self.assertRaisesRegex(ValueError, "intended product"):
            diagnose(case, self.library)

    def test_nominal_with_charge_is_rejected(self):
        case = make_case(51)
        case["observations"][0]["charge"] = 2
        with self.assertRaisesRegex(ValueError, "integer"):
            diagnose(case, self.library)

    def test_exact_trp_from_formula_without_starting_mass(self):
        case = {"reaction": {"description": "global deprotection", "stage_tags": ["acid_global_deprotection", "edt"]},
                "intended_product": {"sequence": "W", "nterm": "free", "cterm": "acid"},
                "observations": [{"kind": "neutral_mono", "value": 376.052704,
                                  "tolerance": 0.000002, "tolerance_unit": "Da"}]}
        result = diagnose(case, self.library)
        self.assertIsNone(result["start_neutral_mono"])
        self.assertEqual("CH03-R009", result["observations"][0]["matches"][0]["events"][0]["rule_id"])

    def test_book_proline_insertion_reconstructs_candidate(self):
        case = {"reaction": {"description": "Fmoc solid phase", "stage_tags": ["fmoc_spps"]},
                "intended_product": {"sequence": "PPPPPPA", "nterm": "free", "cterm": "amide"},
                # H2O + (7 Pro residues) + Ala residue + the terminal amidation delta.
                "observations": [{"kind": "neutral_mono", "value": 767.433010,
                                  "tolerance": 0.000003, "tolerance_unit": "Da"}]}
        events = [e for h in diagnose(case, self.library)["observations"][0]["matches"]
                  for e in h["events"] if e["rule_id"] == "CH10-R005"]
        self.assertEqual({e["site"] for e in events}, {1, 2, 3, 4, 5})
        self.assertTrue(all(e["product_structure"]["sequence"] == list("PPPPPPPA") for e in events))

    def test_book_derived_ions_keep_termini_separate(self):
        for cterm, mz in [("acid", 1059.495294), ("amide", 1058.511279)]:
            with self.subTest(cterm=cterm):
                # Calculated exercise values, not extracted instrument peaks.
                case = {"reaction": {"description": "Fmoc-piperidine", "stage_tags": ["piperidine"]},
                        "intended_product": {"sequence": "DEQEALNSC", "nterm": "free", "cterm": cterm},
                        "observations": [{"kind": "mz", "value": mz, "charge": 1, "ions": {"H": 1},
                                          "tolerance": 0.000003, "tolerance_unit": "Da"}]}
                hit, = diagnose(case, self.library)["observations"][0]["matches"]
                self.assertEqual("CH02-R002", hit["events"][0]["rule_id"])

    def test_book_derived_fmoc_trp_model_ion(self):
        case = {"reaction": {"description": "acid EDT treatment", "stage_tags": ["acid_global_deprotection", "edt"]},
                "intended_product": {"sequence": "W", "nterm": "free", "cterm": "acid",
                                     "modifications": [{"site": "N", "group": "Fmoc"}]},
                "observations": [{"kind": "mz", "value": 599.128060, "charge": 1, "ions": {"H": 1},
                                  "tolerance": 0.000003, "tolerance_unit": "Da"}]}
        hit, = diagnose(case, self.library)["observations"][0]["matches"]
        self.assertEqual("CH03-R009", hit["events"][0]["rule_id"])


if __name__ == "__main__":
    unittest.main()
