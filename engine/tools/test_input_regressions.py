import unittest

from diagnose_case import diagnose, load_rules, structure


def peptide(sequence, **fields):
    return {"sequence": sequence, "nterm": "free", "cterm": "acid", **fields}


class InputRegressions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def result(self, target, mass, *, tags=(), context=None, tolerance=0.000003):
        return diagnose({
            "reaction": {"description": "Regression fixture", "stage_tags": list(tags)},
            "context": context or {},
            "intended_product": target,
            "observations": [{"kind": "neutral_mono", "value": mass,
                              "tolerance": tolerance, "tolerance_unit": "Da"}],
        }, self.library)["observations"][0]

    def event(self, result, rule, bucket="matches"):
        return next(e for hit in result[bucket] for e in hit["events"] if e["rule_id"] == rule)

    def test_trp_unknown_and_known_conditions(self):
        for tags, expected in [((), "unknown"), (("edt",), "unknown"),
                               (("edt", "acid_global_deprotection"), "compatible")]:
            with self.subTest(tags=tags):
                event = self.event(self.result(peptide("W"), 376.052704, tags=tags), "CH03-R009")
                self.assertEqual(event["condition"], expected)

    def test_unrelated_tag_and_wider_tolerance_keep_candidates(self):
        old = self.result(peptide("GC"), 229.142641)["matches"]
        self.assertTrue(old)
        for params in ({"tags": ["oxidation"]}, {"tolerance": 0.001}):
            new = self.result(peptide("GC"), 229.142641, **params)["matches"]
            self.assertTrue(all(hit in new for hit in old))

    def test_explicit_reagent_absence_is_kept_as_conflict(self):
        result = self.result(peptide("GC"), 229.142641, context={"stage_tags_absent": ["piperidine"]})
        self.assertFalse(result["matches"])
        event = self.event(result, "CH02-R002", "condition_conflicts")
        self.assertEqual(event["condition"], "conflict")
        self.assertIn("piperidine", event["conflicting_conditions"])

    def test_contradictory_tags_raise(self):
        with self.assertRaisesRegex(ValueError, "contradictory"):
            self.result(peptide("GC"), 229.142641, tags=["piperidine"],
                        context={"stage_tags_absent": ["piperidine"]})

    def test_equivalent_c_term_amide(self):
        a = peptide("G", cterm="amide")
        b = peptide("G", modifications=[{"site": "C", "group": "C-amide"}])
        self.assertEqual(structure(a)[1], structure(b)[1])
        matches = self.result(a, 75.032028)["matches"]
        self.assertEqual(matches, self.result(b, 75.032028)["matches"])
        self.assertIn("SUP-045", {e["rule_id"] for hit in matches for e in hit["events"]})

    def test_bad_fields_and_incompatible_sites(self):
        for record in (
            peptide("G", modification=[{"site": "N", "group": "Fmoc"}]),
            peptide("G", modifications=[{"site": "N", "group": "Fmoc", "typo": 1}]),
            peptide("G", modifications=[{"site": "N", "group": "Fmoc"}] * 2),
            peptide("G", modifications=[{"site": 1, "group": "Pbf"}]),
            peptide("GK", modifications=[{"site": 2, "group": "Fmoc"}, {"site": 2, "group": "Boc"}]),
        ):
            with self.subTest(record=record), self.assertRaises(ValueError):
                structure(record)
