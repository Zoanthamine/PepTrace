"""Independent elemental and applicability checks for reviewed CH07 branches."""
import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules


class ChapterSevenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def run_case(self, target, observed, tags=(), features=(), absent=()):
        return diagnose({"intended_product": target,
                         "reaction": {"description": "CH07 independent regression",
                                      "stage_tags": list(tags), "stage_tags_absent": list(absent)},
                         "context": {"features_present": list(features)},
                         "observations": [{"kind": "neutral_mono", "value": observed,
                                           "tolerance": 0.000004, "tolerance_unit": "Da"}]},
                        self.library)["observations"][0]

    @staticmethod
    def target(sequence, mods=(), nterm="free", cterm="acid"):
        return dict(sequence=sequence, modifications=list(mods), nterm=nterm, cterm=cterm)

    @staticmethod
    def events(result, rid, bucket="matches"):
        return [e for h in result[bucket] for e in h["events"] if e["rule_id"] == rid]

    def test_transfer_preserves_formula_and_distinct_acceptors(self):
        for rid, sequence, site, group, increment, tags in [
            ("CH07-R002", "GHKK", 2, "Ac", "C2H2O", ["acetyl_transfer"]),
            ("CH07-R010", "SK", 1, "Trifluoroacetyl", "C2H-1F3O", ["base_neutralization"]),
            ("CH07-R012", "GWKK", 2, "Formyl", "CO", ["trp_formyl_deprotection"]),
        ]:
            with self.subTest(rid=rid):
                observed = mass.calculate_mass(sequence=sequence) + mass.calculate_mass(formula=increment)
                target = self.target(sequence, [{"site": site, "group": group}])
                events = self.events(self.run_case(target, observed, tags), rid)
                destinations = {e["site"][1] for e in events}
                self.assertEqual(destinations, {"N"} if rid == "CH07-R010" else {"N", 3, 4})
                self.assertTrue(all(e["delta_formula"] == {} for e in events))
                self.assertTrue(all(e["condition"] == "compatible" for e in events))
                self.assertFalse(self.events(self.run_case(self.target(sequence), observed, tags), rid))

    def test_hydrazide_is_a_distinct_reference(self):
        # G acid C2H5NO2 -> hydrazide C2H7N3O.
        for formula, tag in [("C3H7N3O2", "formic_acid"), ("C4H9N3O2", "acetic_acid")]:
            observed = mass.calculate_mass(formula=formula)
            event, = self.events(self.run_case(self.target("G", cterm="hydrazide"), observed,
                                              [tag, "hydrazide_to_azide"]), "CH07-R014")
            self.assertEqual(event["site"], "C")
            self.assertEqual(event["condition"], "compatible")
            self.assertFalse(self.events(self.run_case(self.target("G"), observed,
                                                       [tag, "hydrazide_to_azide"]), "CH07-R014"))

    def test_nma_cycle_can_be_internal_and_preserves_free_nalpha(self):
        # GKG C10H20N4O4 + Nma C8H7NO + one formaldehyde minus water = +C9H7NO.
        observed = mass.calculate_mass(sequence="GKG") + mass.calculate_mass(formula="C9H7NO")
        target = self.target("GKG", [{"site": 2, "group": "Nma"}])
        event, = self.events(self.run_case(target, observed, ["formaldehyde"]), "CH07-R035")
        self.assertEqual(event["site"], 2)
        self.assertEqual(event["product_structure"]["topology"], "local_cycle")
        self.assertEqual(event["product_structure"]["nterm"], "free")
        self.assertFalse(self.events(self.run_case(self.target("GKG"), observed, ["formaldehyde"]), "CH07-R035"))

    def test_fmoc_loss_selects_terminal_attachment(self):
        target = self.target("GK", [{"site": "N", "group": "Fmoc"}, {"site": 2, "group": "Boc"}])
        observed = mass.calculate_mass(sequence="GK") + mass.calculate_mass(formula="C5H8O2")
        tags = ["alloc_deprotection", "n_silyl_allyl_scavenger"]
        event, = self.events(self.run_case(target, observed, tags), "CH07-R044")
        self.assertEqual(event["site"], "N")
        self.assertEqual(event["product_structure"]["modifications"], [{"site": 2, "group": "Boc"}])
        self.assertFalse(self.events(self.run_case(self.target("GK", [{"site": 2, "group": "Fmoc"}]),
                                                 mass.calculate_mass(sequence="GK"), tags), "CH07-R044"))

    def test_single_bridges_keep_chain_multiplicity_and_atom_sites(self):
        target = self.target("CKR")
        parent = mass.calculate_mass(sequence="CKR")
        tags = ["formaldehyde", "crosslink_favorable_ph"]
        for observed, chains in [(parent + 12, 1), (2 * parent + 12, 2)]:
            events = self.events(self.run_case(target, observed, tags), "CH07-R030")
            cys = [e for e in events if e["site"][0][2] == "Cys_Sgamma"]
            self.assertEqual(len(cys), 3)
            self.assertEqual({e["site"][1][2] for e in cys},
                             {"guanidino_Ndelta", "guanidino_Nomega", "guanidino_Nomega_prime"})
            self.assertTrue(all(len(e["product_structure"]["chains"]) == chains for e in events))
        blocked = self.target("CKR", [{"site": 1, "group": "Acm"}])
        events = self.events(self.run_case(blocked, parent + 12 + mass.calculate_mass(formula="C3H5NO"), tags), "CH07-R030")
        self.assertFalse(any(e["site"][0][2] == "Cys_Sgamma" for e in events))

    def test_external_gly_is_free_molecule_and_isomer_sites_remain(self):
        target = self.target("NQHWY", nterm="acetyl")
        parent = mass.calculate_mass(sequence="NQHWY") + mass.calculate_mass(formula="C2H2O")
        # Free Gly C2H5NO2 + HCHO - H2O = C3H5NO2, 87.032028 Da.
        observed = parent + mass.calculate_mass(formula="C3H5NO2")
        tags = ["formaldehyde", "crosslink_favorable_ph", "free_gly_added"]
        events = self.events(self.run_case(target, observed, tags), "CH07-R030")
        self.assertEqual(len(events), 7)
        self.assertEqual({tuple(e["site"][0][1:]) for e in events},
                         {(1, "sidechain_amide_N"), (2, "sidechain_amide_N"), (3, "imidazole_N"),
                          (4, "indole_N"), (4, "indole_C2"), (5, "phenol_ortho_C3"), (5, "phenol_ortho_C5")})
        self.assertTrue(all(e["delta_formula"] == {"C": 3, "H": 5, "N": 1, "O": 2} for e in events))
        conflicts = self.events(self.run_case(target, observed, tags[:-1], absent=["free_gly_added"]),
                                "CH07-R030", "condition_conflicts")
        self.assertEqual(len(conflicts), 7)
        self.assertFalse(self.events(self.run_case(target, parent + mass.calculate_mass(formula="C2H3NO"), tags), "CH07-R030"))

    def test_allylation_uses_distinct_terminal_and_lys_facts(self):
        observed = mass.calculate_mass(sequence="GK") + mass.calculate_mass(formula="C3H4")
        tags = ["alloc_deprotection", "pd0"]
        events = self.events(self.run_case(self.target("GK"), observed, tags, ["alloc_protected:N"]), "CH07-R040")
        self.assertEqual({e["site"] for e in events if e["condition"] == "compatible"}, {"N"})
        events = self.events(self.run_case(self.target("GK"), observed, tags, ["alloc_protected:2"]), "CH07-R040")
        self.assertEqual({e["site"] for e in events if e["condition"] == "compatible"}, {2})

    def test_book_hypothesis_keeps_evidence_status(self):
        observed = mass.calculate_mass(sequence="GK") + mass.calculate_mass(formula="CO")
        events = self.events(self.run_case(self.target("GK"), observed, ["dmf", "cdi_coupling"]), "CH07-R018")
        self.assertTrue(any(e["source_evidence_level"] == "book_hypothesis" for e in events))
        self.assertTrue(any(e["origin"] == "inference_from_book" for e in events))

    def test_book_trp_iminium_has_fixed_charge_and_neutral_formula_difference(self):
        parent = mass.calculate_mass(sequence="GW")
        # Parent [M+H]+ plus 12 Da; the product replaces indole N-H by N+=CH2.
        observed = parent + 12 + 1.007276466621
        result = diagnose({"intended_product": self.target("GW"),
                           "reaction": {"description": "formaldehyde exposure", "stage_tags": ["formaldehyde"]},
                           "observations": [{"kind": "mz", "value": observed, "charge": 1,
                                             "ions": {}, "intrinsic_charge": 1,
                                             "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)["observations"][0]
        event, = self.events(result, "CH07-R029")
        self.assertEqual(event["delta_formula"], {"C": 1, "H": 1})
        self.assertEqual(event["intrinsic_charge"], 1)
        self.assertEqual(event["site"], 2)
        self.assertFalse(any(e["site"] == 2 for e in self.events(
            self.run_case(self.target("GW"), parent + 12, ["formaldehyde"]), "CH07-R029")))

    def test_beta_ala_monomer_substitution_differs_from_dipeptide_insertion(self):
        # Beta-Ala and alpha-Ala share composition but retain different residue codes.
        parent = mass.calculate_mass(sequence="GA")
        result = self.run_case(self.target("GA"), parent,
                               ["fmoc_osu_beta_ala_monomer_contamination", "coupling"])
        event, = self.events(result, "CH07-R046")
        self.assertEqual(event["site"], 2)
        self.assertEqual(event["product_structure"]["sequence"], ["G", "beta-Ala"])
        self.assertEqual(event["delta_formula"], {})
        self.assertEqual(event["origin"], "inference_from_book")
        self.assertFalse(self.events(result, "CH07-R047"))

    def test_cap_boundaries_preserve_end_groups_and_block_occupied_nterm(self):
        tags = ["acetic_anhydride"]
        parent = mass.calculate_mass(sequence="GK")
        observed = parent + mass.calculate_mass(formula="C2H2O")
        self.assertTrue(self.events(self.run_case(self.target("GK"), observed, tags), "CH07-R001"))
        blocked = self.target("GK", [{"site": "N", "group": "Fmoc"}])
        self.assertFalse(self.events(self.run_case(blocked, observed + mass.calculate_mass(formula="C15H10O2"), tags), "CH07-R001"))
        # Mid-synthesis cap uses the retained K amide rather than the parent GK amide.
        observed = mass.calculate_mass(formula="C8H17N3O2")
        event, = self.events(self.run_case(self.target("GK", cterm="amide"), observed, tags,
                                          ["synthesis_boundary:2"]), "CH07-R001")
        self.assertEqual(event["product_structure"]["sequence"], ["K"])
        self.assertEqual(event["condition"], "compatible")


if __name__ == "__main__":
    unittest.main()
