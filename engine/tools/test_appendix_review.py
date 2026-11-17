import json
import unittest

from pyteomics import mass

from diagnose_case import BOOK, diagnose


class AppendixReviewTests(unittest.TestCase):
    def candidate(self, sequence, delta, tags=(), context=None):
        return diagnose({
            "intended_product": {"sequence": sequence, "nterm": "free", "cterm": "acid"},
            "reaction": {"description": "appendix scheme", "stage_tags": list(tags)},
            "context": context or {},
            "observations": [{"kind": "neutral_mono", "value": mass.calculate_mass(sequence=sequence) + delta,
                              "tolerance": 0.000003, "tolerance_unit": "Da"}],
        })["observations"][0]

    def test_cyanation_isobars_have_distinct_sites_and_connections(self):
        result = self.candidate("GMCG", 24.995249, ["thiol_cyanation", "iodoacetonitrile"])
        events = {e["rule_id"]: e for h in result["matches"] for e in h["events"]}
        for rule, site, group in [("APP-I-R035", 2, "MetCyanoMethylExchange"), ("APP-I-R036", 3, "S-Cyano")]:
            event = events[rule]
            self.assertEqual(event["site"], site)
            self.assertEqual(event["origin"], "book_appendix")
            self.assertEqual(event["product_structure"]["modifications"][0]["group"], group)
            self.assertNotIn("intrinsic_charge", event)

    def test_met_degradation_matches_figure_not_literal_name(self):
        hits = self.candidate("GMG", -48.003371, ["haloacetyl_induced_met_degradation"])["matches"]
        event, = [e for h in hits for e in h["events"] if e["rule_id"] == "APP-I-R004"]
        self.assertEqual(event["product_structure"]["sequence"], ["G", "VinylGly", "G"])
        self.assertIn("allylglycine", event["diagnostic_evidence"][0])

    def test_appendix_local_ring_requirements(self):
        cases = [
            ("APP-I-R039", "AGA", "APA", 25.979265, [], ["hydantoin_carbonyl_source"]),
            ("APP-I-R052", "ASA", "SAA", 40.031300, ["global_deprotection"], ["pseudoproline_used:2"]),
        ]
        for rid, positive, negative, delta, tags, features in cases:
            with self.subTest(rid=rid):
                for sequence, present in [(positive, True), (negative, False)]:
                    result = self.candidate(sequence, delta, tags, {"features_present": features})
                    self.assertEqual(any(e["rule_id"] == rid and e["condition"] == "compatible"
                                         for h in result["matches"] for e in h["events"]), present)

    def test_fixed_charge_ledger_excludes_neutral_mass(self):
        rows = json.loads((BOOK / "appendix_I_element_account.json").read_text(encoding="utf-8"))["rows"]
        row = next(r for r in rows if r["id"] == "AI-072")
        self.assertIsNone(row["exact_neutral_delta"])
        self.assertEqual(row["intrinsic_charge"], 1)
        self.assertAlmostEqual(row["elemental_mass_delta"], mass.calculate_mass(formula="C6H13S"), places=5)
        self.assertAlmostEqual(row["same_protonated_ion_delta"], mass.calculate_mass(formula="C6H12S"), places=5)
