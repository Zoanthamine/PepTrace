import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules, structure


class ChapterNineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    @staticmethod
    def target(seq, mods=(), bonds=()):
        value = dict(sequence=seq, nterm="free", cterm="acid", modifications=list(mods))
        if bonds:
            value.update(bonds=list(bonds), topology="crosslinked")
        return value

    def run_case(self, target, value, tags=(), absent=()):
        return diagnose({"intended_product": target,
                         "reaction": {"description": "CH09 independent regression", "stage_tags": list(tags), "stage_tags_absent": list(absent)},
                         "observations": [{"kind": "neutral_mono", "value": value, "tolerance": 0.000004, "tolerance_unit": "Da"}]}, self.library)

    @staticmethod
    def events(result, n, bucket="matches"):
        return [e for h in result["observations"][0][bucket] for e in h["events"] if e["rule_id"] == f"CH09-R{n:03}"]

    def test_all_unmodified_local_products_have_independent_elemental_deltas(self):
        for n, residue, formula, count in [
            (1, "C", "O", 1), (2, "C", "O2", 1), (3, "C", "O3", 1),
            (29, "M", "O", 2), (30, "M", "O2", 1), (31, "M", "O", 1),
            (36, "W", "O", 1), (37, "W", "O2", 1), (38, "W", "O3", 1),
            (39, "W", "C-1O", 1), (40, "W", "C-1O2", 1), (41, "W", "O", 1), (42, "W", "O", 1),
            (44, "H", "O", 1), (45, "Y", "O", 1), (46, "Y", "H-1NO2", 1),
            (47, "Y", "H-1Cl", 1), (48, "Y", "H-2Cl2", 1), (49, "Y", "H-1Br", 1), (50, "Y", "H-2Br2", 1),
            (52, "F", "O", 1), (53, "F", "O", 1), (54, "F", "O2", 1), (57, "W", "H2", 1),
        ]:
            with self.subTest(rule=n):
                seq = "G"+residue+"G"
                observed = mass.calculate_mass(sequence=seq)+mass.calculate_mass(formula=formula)
                found = self.events(self.run_case(self.target(seq), observed), n)
                self.assertEqual(len(found), count)
                self.assertEqual({e["site"] for e in found}, {2})
                self.assertFalse(self.events(self.run_case(self.target("GAG"), mass.calculate_mass(sequence="GAG")+mass.calculate_mass(formula=formula)), n))

    def test_disulfide_inputs_and_oxidation_states_are_explicit(self):
        free = mass.calculate_mass(sequence="GCCG")
        hydrogen = mass.calculate_mass(formula="H2")
        oxygen = mass.calculate_mass(formula="O")
        ss = {"kind": "disulfide", "sites": [2, 3]}
        target = self.target("GCCG", bonds=[ss])
        self.assertAlmostEqual(mass.calculate_mass(composition=structure(target)[1]), free-hydrogen, places=7)
        result = self.run_case(target, free-hydrogen+oxygen)
        found = self.events(result, 6)
        self.assertEqual(len(found), 2)
        self.assertEqual({e["bond_edit"]["reacting_site"] for e in found}, {2, 3})
        self.assertTrue(result["scope"]["sequence_rules_skipped"])
        full_oxidation = self.events(self.run_case(target, free-hydrogen+2*oxygen), 7)
        self.assertEqual(len(full_oxidation), 2)
        self.assertTrue(all(e["delta_formula"] == {"O": 2} for e in full_oxidation))
        for e in found:
            oxidized = e["product_structure"]
            # Produced structures can be validated and used as the next explicit precursor.
            self.assertAlmostEqual(mass.calculate_mass(composition=structure(oxidized)[1]), free-hydrogen+oxygen, places=7)
            event, = self.events(self.run_case(oxidized, free-hydrogen+2*oxygen), 7)
            self.assertEqual(event["delta_formula"], {"O": 1})
        with self.assertRaisesRegex(ValueError, "oxidized_site"):
            structure(self.target("GCCG", bonds=[{"kind": "thiosulfinate", "sites": [2, 3]}]))
        with self.assertRaisesRegex(ValueError, "overlaps"):
            structure(self.target("GCCG", [{"site": 2, "group": "Acm"}], [ss]))

    def test_bridge_hydrolysis_oxidative_cleavage_and_reduction_mass_balance(self):
        free = mass.calculate_mass(sequence="GCCG")
        ss = self.target("GCCG", bonds=[{"kind": "disulfide", "sites": [2, 3]}])
        tso = self.target("GCCG", bonds=[{"kind": "thiosulfinate", "sites": [2, 3], "oxidized_site": 2}])
        event, = self.events(self.run_case(tso, free+mass.calculate_mass(formula="O2")), 9)
        self.assertEqual(event["delta_formula"], {"H": 2, "O": 1})
        self.assertEqual(event["product_structure"]["modifications"], [{"site": 2, "group": "CysSulfenic"}, {"site": 3, "group": "CysSulfenic"}])
        event, = self.events(self.run_case(ss, free+mass.calculate_mass(formula="O6")), 10)
        self.assertEqual(event["delta_formula"], {"H": 2, "O": 6})
        event, = self.events(self.run_case(ss, free, ["dtt"]), 58)
        self.assertEqual(event["delta_formula"], {"H": 2})
        self.assertNotIn("bonds", event["product_structure"])

    def test_disulfide_beta_elimination_preserves_whole_chain_composition(self):
        target = self.target("GCCG", bonds=[{"kind": "disulfide", "sites": [2, 3]}])
        observed = mass.calculate_mass(sequence="GCCG")-mass.calculate_mass(formula="H2")
        events = self.events(self.run_case(target, observed, ["base_exposure"]), 25)
        self.assertEqual(len(events), 2)
        self.assertEqual({tuple(e["product_structure"]["sequence"]) for e in events},
                         {("G", "Dha", "C", "G"), ("G", "C", "Dha", "G")})
        self.assertTrue(all(e["delta_formula"] == {} for e in events))
        self.assertTrue(all("bonds" not in e["product_structure"] for e in events))

    def test_intrachain_and_interchain_disulfide_or_dityrosine_masses(self):
        free = mass.calculate_mass(sequence="GCCG")
        loss = mass.calculate_mass(formula="H2")
        event, = self.events(self.run_case(self.target("GCCG"), free-loss), 4)
        self.assertEqual(len(event["product_structure"]["chains"]), 1)
        found = self.events(self.run_case(self.target("GCCG"), 2*free-loss), 5)
        self.assertEqual(len(found), 3)
        self.assertTrue(all(len(e["product_structure"]["chains"]) == 2 for e in found))
        parent = mass.calculate_mass(sequence="GYYG")
        for observed, count in [(parent-loss, 1), (2*parent-loss, 3)]:
            self.assertEqual(len(self.events(self.run_case(self.target("GYYG"), observed), 51)), count)

    def test_precursor_groups_and_sulfenyl_amide_backbone_hydrogen(self):
        parent = mass.calculate_mass(sequence="GCG")
        for n, group, observed, delta in [
            (18, "CysSulfenic", parent, {"O": -1}), (21, "CysSulfinic", parent, {"O": -2}),
            (24, "CysSulfonic", parent-mass.calculate_mass(formula="H2S"), {"H": -2, "O": -3, "S": -1}),
            (26, "Persulfide", parent+mass.calculate_mass(formula="SO3"), {"O": 3}),
            (15, "CysSulfenic", parent-mass.calculate_mass(formula="H2"), {"H": -2, "O": -1}),
            (16, "CysSulfenic", parent+mass.calculate_mass(formula="C6HN3O4"), {"C": 6, "H": 1, "N": 3, "O": 3}),
        ]:
            with self.subTest(rule=n):
                event, = self.events(self.run_case(self.target("GCG", [{"site": 2, "group": group}]), observed), n)
                self.assertEqual(event["delta_formula"], delta)
                self.assertFalse(self.events(self.run_case(self.target("GCG"), observed), n))
        self.assertFalse(self.events(self.run_case(self.target("GCPG", [{"site": 2, "group": "CysSulfenic"}]),
            mass.calculate_mass(sequence="GCPG")-mass.calculate_mass(formula="H2")), 15))

    def test_protected_thioether_oxidation_retains_correct_group(self):
        parent = mass.calculate_mass(sequence="GCG")
        for group, formula in [("Acm", "C3H5NO"), ("Bzl", "C7H6")]:
            target = self.target("GCG", [{"site": 2, "group": group}])
            event, = self.events(self.run_case(target, parent+mass.calculate_mass(formula=formula)+mass.calculate_mass(formula="O")), 27)
            self.assertEqual(event["delta_formula"], {"O": 1})
            self.assertEqual(len(event["product_structure"]["modifications"]), 1)

    def test_met_reduction_associated_disulfide_requires_met_sequence(self):
        free = mass.calculate_mass(sequence="MCCG")
        event, = self.events(self.run_case(self.target("MCCG"), free-mass.calculate_mass(formula="H2")), 33)
        self.assertEqual(event["site"], [[1, 2, "Cys_Sgamma"], [1, 3, "Cys_Sgamma"]])
        self.assertFalse(self.events(self.run_case(self.target("GCCG"), mass.calculate_mass(sequence="GCCG")-mass.calculate_mass(formula="H2")), 33))

    def test_trp_dimerization_and_nfk_reference(self):
        parent = mass.calculate_mass(sequence="GWG")
        for observed in [2*parent, 2*parent-mass.calculate_mass(formula="H2")]:
            event, = self.events(self.run_case(self.target("GWG"), observed), 34)
            self.assertEqual(len(event["product_structure"]["chains"]), 2)
        nfk = self.target(["G", "NFK", "G"])
        observed = mass.calculate_mass(formula="C14H18N4O5")  # two Gly residues + Kyn residue + H2O
        event, = self.events(self.run_case(nfk, observed, ["deformylation"]), 39)
        self.assertEqual(event["delta_formula"], {"C": -1, "O": -1})

    def test_explicit_azido_precursors_and_conflicting_literature(self):
        for number, precursor, product in [(59, "Nva(N3)", "Orn"), (60, "Nle(N3)", "K")]:
            target = self.target(["G", precursor, "G"])
            formula = "C9H18N4O4" if product == "Orn" else "C10H20N4O4"
            event, = self.events(self.run_case(target, mass.calculate_mass(formula=formula)), number)
            self.assertEqual(event["delta_formula"], {"H": 2, "N": -2})
            if number == 60:
                self.assertEqual(event["source_evidence_level"], "book_conflicting_reports")
        target = self.target("GKG", [{"site": "N", "group": "NalphaAzide"}])
        event, = self.events(self.run_case(target, mass.calculate_mass(sequence="GKG")), 61)
        self.assertEqual(event["site"], "N")
        self.assertFalse(self.events(self.run_case(self.target("GKG"), mass.calculate_mass(sequence="GKG")), 61))
        with self.assertRaisesRegex(ValueError, "primary N-terminal amine"):
            structure(self.target("PG", [{"site": "N", "group": "NalphaAzide"}]))


if __name__ == "__main__":
    unittest.main()
