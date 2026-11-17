import unittest

from pyteomics import mass

from diagnose_case import diagnose, load_rules


class LinearRuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    def results(self, sequence, value, *, cterm="acid", nterm="free", tags=(), context=None):
        case = {
            "intended_product": {"sequence": sequence, "nterm": nterm, "cterm": cterm},
            "reaction": {"description": "Chapter 1 structural test", "stage_tags": list(tags)},
            "context": context or {},
            "observations": [{"kind": "neutral_mono", "value": value,
                              "tolerance": 0.000003, "tolerance_unit": "Da"}],
        }
        return diagnose(case, self.library)["observations"][0]

    def events(self, rule, *args, bucket="matches", **kwargs):
        return [event for hit in self.results(*args, **kwargs)[bucket]
                for event in hit["events"] if event["rule_id"] == rule]

    def test_sequence_rule_exclusion_is_visible(self):
        case = {"intended_product": {"sequence": "GDA", "nterm": "free", "cterm": "acid",
                "modifications": [{"site": 1, "delta_formula": {"C": 1}}]},
                "reaction": {"description": "unknown modification"},
                "observations": [{"kind": "neutral_mono", "value": 300,
                                  "tolerance": 1, "tolerance_unit": "Da"}]}
        scope = diagnose(case, self.library)["scope"]
        self.assertTrue(scope["sequence_rules_skipped"])
        self.assertEqual(scope["sequence_rules_blocking_modifications"][0]["site"], 1)
        self.assertTrue(scope["sequence_rules_skip_reason"])

    def test_ch07_side_chain_acylation_mass_and_sites(self):
        sequence = "GKSR"
        parent = mass.calculate_mass(sequence=sequence)
        acetyl = mass.calculate_mass(formula="C2H2O")
        lys, = self.events("CH07-R003", sequence, parent + acetyl, tags=["acetic_anhydride"])
        arg, = self.events("CH07-R004", sequence, parent + acetyl, tags=["acetic_anhydride"])
        self.assertEqual(lys["site"], 2)
        self.assertEqual(arg["site"], 4)
        tfa = mass.calculate_mass(formula="C2F3O") - mass.calculate_mass(formula="H")
        ser, = self.events("CH07-R009", sequence, parent + tfa, tags=["tfa"])
        self.assertEqual(ser["site"], 3)
        self.assertFalse(self.events("CH07-R009", "GKA", mass.calculate_mass(sequence="GKA") + tfa))

    def test_ch07_trifluoroacetylated_b_segment_is_not_a_full_length_adduct(self):
        # AG|KST ligation, intended C-terminal amide. Formula from KST amide plus the acyl group.
        # Free KST: C13H26N4O6; acid->amide: +NH-O (H increases by one);
        # trifluoroacetylation: +C2F3O-H. Independent elemental ledger.
        suffix = mass.calculate_mass(formula="C15H26F3N5O6")
        found, = self.events("CH07-R006", "AGKST", suffix, cterm="amide",
                             tags=["tfa", "segment_condensation"],
                             context={"features_present": ["segment_b_start:3"]})
        self.assertEqual(found["product_structure"]["sequence"], list("KST"))
        self.assertEqual(found["product_structure"]["cterm"], "amide")
        self.assertEqual(found["condition"], "compatible")
        excluded, = self.events("CH07-R006", "AGKST", suffix, cterm="amide",
                                tags=["tfa", "segment_condensation"],
                                context={"features_absent": ["segment_b_start:3"]},
                                bucket="condition_conflicts")
        self.assertEqual(excluded["conflicting_conditions"], ["segment_b_start:3"])

    def test_ch07_site_specific_transform_and_isobars(self):
        parent = mass.calculate_mass(sequence="KGS")
        tfa = 95.982298
        positions = {e["site"] for e in self.events("CH07-R006", "KGS", parent + tfa, tags=["tfa"])}
        self.assertEqual(positions, {"N", 1})
        events = self.events("CH07-R006", "KGS", parent + tfa + 42.010565, tags=["tfa"], nterm="acetyl")
        self.assertEqual([e["site"] for e in events], [1])
        formyl = mass.calculate_mass(formula="CO")
        found = self.events("CH07-R013", "WKST", mass.calculate_mass(sequence="WKST") + formyl,
                            tags=["formic_acid"])
        self.assertEqual({e["site"] for e in found}, {"N", 1, 2, 3, 4})
        propyl = mass.calculate_mass(sequence="PW") + 12
        self.assertFalse(self.events("CH07-R029", "PW", propyl, tags=["formaldehyde"]))
        self.assertFalse(self.events("CH07-R032", "GP", mass.calculate_mass(sequence="GP") + 12))

    def test_ch07_cbz_conversion_depends_on_reference_and_preserves_other_sites(self):
        target = {"sequence": "KGS", "nterm": "free", "cterm": "acid",
                  "modifications": [{"site": 1, "group": "Cbz"}]}
        expected = mass.calculate_mass(sequence="KGS") + mass.calculate_mass(formula="C7H6")
        result = diagnose({"intended_product": target,
                           "reaction": {"description": "Boc cleavage", "stage_tags": ["tfa", "boc_deprotection"]},
                           "observations": [{"kind": "neutral_mono", "value": expected,
                                             "tolerance": 0.000003, "tolerance_unit": "Da"}]}, self.library)
        event, = [e for h in result["observations"][0]["matches"] for e in h["events"] if e["rule_id"] == "CH07-R026"]
        self.assertEqual(event["delta_formula"], {"C": -1, "O": -2})
        self.assertEqual(event["site"], 1)
        self.assertFalse(self.events("CH07-R026", "KGS", expected, tags=["tfa", "boc_deprotection"]))

    def test_asp_pro_independent_fragments_from_two_sites(self):
        sequence = "GDPDPK"
        for fragment, site, c_amide in [("GD", [2, 3], False), ("PDPK", [2, 3], True),
                                       ("GDPD", [4, 5], False), ("PK", [4, 5], True)]:
            with self.subTest(fragment=fragment):
                expected = mass.calculate_mass(sequence=fragment)
                if c_amide:
                    expected += mass.calculate_mass(formula="NH3") - mass.calculate_mass(formula="H2O")
                event, = self.events("CH01-R008", sequence, expected, cterm="amide")
                self.assertEqual(event["site"], site)
                self.assertEqual(event["product_structure"]["sequence"], list(fragment))
                self.assertEqual(event["product_structure"]["cterm"], "amide" if c_amide else "acid")
        self.assertFalse(self.events("CH01-R008", "GNPDPK", mass.calculate_mass(sequence="GN")))

    def test_isodipeptide_loss_same_mass_different_sites(self):
        predicted = mass.calculate_mass(sequence="GSA")
        events = self.events("CH01-R002", "GSSA", predicted, tags=["o_acyl_isodipeptide_activation"],
                             context={"features_present": ["o_acyl_isodipeptide_site:2"],
                                      "features_absent": ["o_acyl_isodipeptide_site:3"]})
        self.assertEqual([(e["site"], e["condition"]) for e in events], [(2, "compatible")])
        conflicts = self.events("CH01-R002", "GSSA", predicted, tags=["o_acyl_isodipeptide_activation"],
                                context={"features_absent": ["o_acyl_isodipeptide_site:3"]},
                                bucket="condition_conflicts")
        self.assertEqual([e["site"] for e in conflicts], [3])
        self.assertEqual({e["site"] for e in self.events("CH01-R002", "GSSA", predicted)}, {2, 3})

    def test_dkp_short_peptide_and_cyclic_coproduct(self):
        for rule in ("CH01-R009", "CH01-R018"):
            short, = self.events(rule, "HPG", mass.calculate_mass(sequence="G"))
            self.assertEqual(short["product_structure"]["sequence"], ["G"])
            ring, = self.events(rule, "HPG", mass.calculate_mass(formula="C11H14N4O2"))
            self.assertEqual(ring["product_structure"]["topology"], "head_tail_cycle")
            self.assertFalse(self.events(rule, "HPG", mass.calculate_mass(sequence="G"), nterm="acetyl"))

    def test_orn_substitution_requires_arg(self):
        predicted = mass.calculate_mass(formula="C7H15N3O3")
        event, = self.events("CH01-R015", "GR", predicted,
                             tags=["guanidino_acylation_risk", "base_after_coupling"])
        self.assertEqual(event["product_structure"]["sequence"], ["G", "Orn"])
        self.assertEqual(event["delta_formula"], {"C": -1, "H": -2, "N": -2})
        self.assertEqual(event["condition"], "compatible")
        self.assertFalse(self.events("CH01-R015", "GK", predicted))

    def test_phosphate_branches_use_phosphorylated_target(self):
        target = {"sequence": "GST", "nterm": "free", "cterm": "acid",
                  "modifications": [{"site": 2, "group": "Phospho"}, {"site": 3, "group": "Phospho"}]}
        parent = mass.calculate_mass(sequence="GST") + 2 * mass.calculate_mass(formula="HPO3")
        loss = mass.calculate_mass(formula="H3PO4")
        for rule, amine, tag, expected_sites in [
            ("CH02-R006", None, "piperidine", {2, 3}),
            ("CH02-R007", "C5H11N", "piperidine", {2, 3}),
            ("CH02-R008", "C6H13N", "4-methylpiperidine", {2, 3}),
            ("CH02-R009", "C4H9N", "pyrrolidine", {2}),
        ]:
            with self.subTest(rule=rule):
                observed = parent - loss + (mass.calculate_mass(formula=amine) if amine else 0)
                case = {"intended_product": target, "reaction": {"description": "phosphate", "stage_tags": [tag]},
                        "observations": [{"kind": "neutral_mono", "value": observed,
                                          "tolerance": 0.000003, "tolerance_unit": "Da"}]}
                output = diagnose(case, self.library)["observations"][0]
                events = [e for h in output["matches"] for e in h["events"] if e["rule_id"] == rule]
                self.assertEqual({e["site"] for e in events}, expected_sites)
                for e in events:
                    self.assertNotIn(e["site"], {m["site"] for m in e["product_structure"]["modifications"]})
                unmodified = dict(case, intended_product=dict(target, modifications=[]))
                self.assertFalse(any(e["rule_id"] == rule for h in
                    diagnose(unmodified, self.library)["observations"][0]["matches"] for e in h["events"]))

    def test_phosphine_elimination_requires_stbu_history_for_support(self):
        expected = mass.calculate_mass(sequence="GC") - mass.calculate_mass(formula="H2S")
        event, = self.events("CH02-R004", "GC", expected, tags=["tributylphosphine"])
        self.assertEqual(event["condition"], "unknown")
        confirmed, = self.events("CH02-R004", "GC", expected, tags=["tributylphosphine"],
                                 context={"features_present": ["cys_stbu_at_reaction_site:2"]})
        self.assertEqual(confirmed["condition"], "compatible")

    def test_absolute_and_signed_relative_observations_match_same_products(self):
        target = {"sequence": "GSSA", "nterm": "free", "cterm": "acid"}
        case = {"intended_product": target, "starting_material": dict(target, nterm="acetyl"),
                "reaction": {"description": "Independent mass reference comparison"}}
        product = mass.calculate_mass(sequence="GSA")
        target_mass = mass.calculate_mass(sequence="GSSA")
        acetyl = mass.calculate_mass(formula="C2H2O")
        observations = [
            {"kind": "neutral_mono", "value": product, "tolerance": 0.000003, "tolerance_unit": "Da"},
            {"kind": "delta_mono", "value": product - target_mass, "reference": "intended_product",
             "tolerance": 0.000003, "tolerance_unit": "Da"},
            {"kind": "delta_mono", "value": product - target_mass - acetyl, "reference": "starting_material",
             "tolerance": 0.000003, "tolerance_unit": "Da"},
            {"kind": "delta_approx", "value": -87, "reference": "intended_product",
             "tolerance": 0.2, "tolerance_unit": "Da"},
        ]
        for observation in observations:
            with self.subTest(observation=observation):
                result = diagnose(dict(case, observations=[observation]), self.library)["observations"][0]
                hits = [hit for hit in result["matches"] if hit["events"][0]["rule_id"] == "CH01-R002"]
                self.assertEqual({hit["events"][0]["site"] for hit in hits}, {2, 3})
                self.assertTrue(all(abs(hit["neutral_mono"] - product) < 0.000003 for hit in hits))
                if observation["kind"] != "neutral_mono":
                    self.assertTrue(all(hit["err_ppm"] is None for hit in hits))

    def test_trp_tert_butyl_and_n_carboxy_create_distinct_products(self):
        parent = mass.calculate_mass(sequence="GWW")
        for rule, net, feature, sites in [
            ("CH03-R001", "C4H8", "tert_butyl_electrophile_available", {2, 3}),
            ("CH03-R002", "CO2", "trp_indole_boc_at_reaction_site:2", {2, 3}),
        ]:
            with self.subTest(rule=rule):
                predicted = parent + mass.calculate_mass(formula=net)
                events = self.events(rule, "GWW", predicted, tags=["acid_global_deprotection"],
                                     context={"features_present": [feature]})
                self.assertEqual({e["site"] for e in events}, sites)
                self.assertTrue(all(e["condition"] == ("compatible" if rule == "CH03-R001" or e["site"] == 2
                                                       else "unknown") for e in events))
                self.assertFalse(self.events(rule, "GFF", mass.calculate_mass(sequence="GFF") +
                                             mass.calculate_mass(formula=net)))

    def test_linker_additions_and_c_terminal_isobar(self):
        amide = mass.calculate_mass(formula="NH3") - mass.calculate_mass(formula="H2O")
        for rule, formula, resin in [("CH03-R003", "C7H6O", "rink_resin"),
                                     ("CH03-R004", "C14H19NO4", "pal_mbha_resin"),
                                     ("CH03-R005", "C9H9NO2", "rink_mbha_resin"),
                                     ("CH03-R006", "C9H5F3O2", "wang_resin")]:
            with self.subTest(rule=rule):
                expected = mass.calculate_mass(sequence="GW") + amide + mass.calculate_mass(formula=formula)
                tags = ["acid_global_deprotection", "tfa", resin]
                hits = self.results("GW", expected, cterm="amide", tags=tags)["matches"]
                reported = {(e["rule_id"], e["site"]) for hit in hits for e in hit["events"]}
                self.assertIn((rule, 2), reported)
                if rule == "CH03-R003":
                    self.assertIn(("CH03-R008", "C"), reported)
        expected = mass.calculate_mass(sequence="GK") + amide + mass.calculate_mass(formula="C7H6O")
        event, = self.events("CH03-R008", "GK", expected, cterm="amide",
                             tags=["acid_global_deprotection", "rink_resin"])
        self.assertEqual(event["product_structure"]["cterm"], "n_4_hydroxybenzyl_amide")
        self.assertFalse(self.events("CH03-R008", "GK", expected))

    def test_bpa_edt_uses_ketone_substrate(self):
        expected = mass.calculate_mass(formula="C20H22N2O3S2")
        # Gly-Bpa = C18H18N2O4. Addition of EDT minus water gives C20H22N2O3S2.
        event, = self.events("CH03-R010", ["G", "Bpa"], expected, tags=["edt", "acid_global_deprotection"])
        self.assertEqual(event["delta_formula"], {"C": 2, "H": 4, "S": 2, "O": -1})
        self.assertFalse(self.events("CH03-R010", "GF", expected))

    def test_trp_crosslink_multiplicity_and_sites(self):
        parent = mass.calculate_mass(sequence="GWW")
        hydrogen = mass.calculate_mass(formula="H2")
        reduced = self.events("CH03-R011", "GWW", 2 * parent, tags=["acid_global_deprotection"])
        self.assertEqual(len(reduced), 4)
        oxidized = self.events("CH03-R012", "GWW", 2 * parent - hydrogen,
                               tags=["acid_global_deprotection", "ddq"])
        self.assertEqual(len(oxidized), 3)
        self.assertTrue(all(len(event["product_structure"]["chains"]) == 2 for event in oxidized))
        intra, = self.events("CH03-R012", "GWW", parent - hydrogen,
                             tags=["acid_global_deprotection", "ddq"])
        self.assertEqual(intra["site"], [[1, 2], [1, 3]])
        self.assertEqual(intra["delta_formula"], {"H": -2})
        self.assertFalse(self.events("CH03-R012", "GW", mass.calculate_mass(sequence="GW") - hydrogen))

    def test_dodt_sulfonium_is_separate_from_neutral_mass(self):
        case = {
            "intended_product": {"sequence": "GM", "nterm": "free", "cterm": "acid"},
            "reaction": {"description": "DODT cleavage", "stage_tags": ["acid_global_deprotection", "dodt"]},
            "context": {"features_present": ["tert_butyl_electrophile_available"]},
        }
        # GM is C7H14N2O3S. The bare sulfonium is C13H27N2O3S2+, electron removed.
        atomic_mass = mass.calculate_mass(formula="C13H27N2O3S2")
        inputs = [
            {"kind": "mz", "value": atomic_mass - 0.000548579909, "charge": 1,
             "intrinsic_charge": 1, "ions": {}},
            {"kind": "mz", "value": (atomic_mass - 0.000548579909 + 1.007276466621) / 2, "charge": 2,
             "intrinsic_charge": 1, "ions": {"H": 1}},
        ]
        for observation in inputs:
            observation.update(tolerance=0.000003, tolerance_unit="Da")
            hit, = diagnose(dict(case, observations=[observation]), self.library)["observations"][0]["matches"]
            self.assertEqual(hit["events"][0]["rule_id"], "CH03-R024")
            self.assertIsNone(hit["neutral_mono"])
            self.assertIsNone(hit["delta_vs_target"])
        neutral = self.events("CH03-R024", "GM", atomic_mass)
        self.assertFalse(neutral)

    def test_acm_loss_and_isomeric_migration(self):
        parent = mass.calculate_mass(sequence="GYC") + mass.calculate_mass(formula="C3H5NO")
        case = {"intended_product": {"sequence": "GYC", "nterm": "free", "cterm": "acid",
                                     "modifications": [{"site": 3, "group": "Acm"}]},
                "reaction": {"description": "Acm intended to remain", "stage_tags": ["acid_global_deprotection"]},
                "observations": [{"kind": "neutral_mono", "value": parent, "tolerance": 0.000003, "tolerance_unit": "Da"},
                                 {"kind": "neutral_mono", "value": mass.calculate_mass(sequence="GYC"),
                                  "tolerance": 0.000003, "tolerance_unit": "Da"}]}
        isomers, loss = diagnose(case, self.library)["observations"]
        self.assertTrue(any(hit["kind"] == "intended_product" for hit in isomers["matches"]))
        migration, = [event for hit in isomers["matches"] for event in hit["events"]
                      if event["rule_id"] == "CH03-R023"]
        self.assertEqual(migration["site"], [2, 3])
        self.assertEqual(migration["delta_formula"], {})
        self.assertTrue(any(e["rule_id"] == "CH03-R021" for hit in loss["matches"] for e in hit["events"]))

    def test_backbone_ester_isomer_and_hydrolysis(self):
        parent = mass.calculate_mass(sequence="GSK")
        rearranged, = self.events("CH04-R001", "GSK", parent, tags=["acid_exposure"])
        self.assertEqual(rearranged["product_structure"]["topology"], "depsipeptide")
        self.assertEqual(rearranged["site"], [1, 2])
        self.assertEqual(rearranged["delta_formula"], {})
        self.assertFalse(self.events("CH04-R001", "SGK", parent))
        for fragment in ("G", "SK"):
            event, = self.events("CH04-R002", "GSK", mass.calculate_mass(sequence=fragment),
                                 tags=["acid_exposure", "aqueous_exposure"])
            self.assertEqual(event["product_structure"]["sequence"], list(fragment))
            self.assertEqual(event["condition"], "compatible")

    def test_terminal_acyl_migration_loss_and_trifluoroacetylation(self):
        parent = mass.calculate_mass(sequence="SG")
        acetyl = mass.calculate_mass(formula="C2H2O")
        shifted, = self.events("CH04-R003", "SG", parent + acetyl, nterm="acetyl")
        self.assertEqual(shifted["product_structure"]["modifications"][0]["site"], 1)
        removed, = self.events("CH04-R004", "SG", parent, nterm="acetyl")
        self.assertEqual(removed["product_structure"]["nterm"], "free")
        n_tfa = mass.calculate_mass(formula="C2F3O") - mass.calculate_mass(formula="H")
        product, = self.events("CH04-R005", "SG", parent + n_tfa, tags=["tfa", "neutralization"])
        self.assertEqual(product["site"], "N")
        self.assertFalse(self.events("CH04-R005", "GS", parent + n_tfa))

    def test_ser_repeated_insertion_and_his_gly_insertion(self):
        predicted = mass.calculate_mass(sequence="GSSK")
        event, = self.events("CH04-R008", "GSK", predicted, tags=["minimal_protection", "base_deprotection"])
        self.assertEqual(event["product_structure"]["sequence"], list("GSSK"))
        self.assertEqual(event["condition"], "compatible")
        expected = mass.calculate_mass(sequence="GAGHK")
        events = self.events("CH04-R010", "GAHK", expected,
                             tags=["his_imidazole_unprotected", "gly_coupling", "base_deprotection"])
        self.assertEqual(len(events), 3)
        self.assertTrue(all("G" in e["product_structure"]["sequence"] for e in events))

    def test_phosphorylated_repeated_ser(self):
        target = {"sequence": "GSK", "nterm": "free", "cterm": "acid",
                  "modifications": [{"site": 2, "group": "Phospho"}]}
        expected = mass.calculate_mass(sequence="GSSK") + 2 * mass.calculate_mass(formula="HPO3")
        case = {"intended_product": target, "reaction": {"description": "global phosphorylation",
                "stage_tags": ["minimal_protection", "global_phosphorylation"]},
                "observations": [{"kind": "neutral_mono", "value": expected,
                                  "tolerance": 0.000003, "tolerance_unit": "Da"}]}
        result = diagnose(case, self.library)["observations"][0]
        event, = [e for h in result["matches"] for e in h["events"] if e["rule_id"] == "CH04-R009"]
        self.assertEqual({m["site"] for m in event["product_structure"]["modifications"]}, {2, 3})

    def test_activation_cyanation_and_residue_insertions(self):
        dehydration = mass.calculate_mass(formula="H2O")
        for rule, sequence, tag in [("CH05-R011", "GNK", "asn_carboxyl_activation"),
                                    ("CH05-R012", "GQK", "gln_carboxyl_activation")]:
            observed = mass.calculate_mass(sequence=sequence) - dehydration
            event, = self.events(rule, sequence, observed, tags=[tag])
            self.assertEqual(event["site"], 2)
            self.assertEqual(event["delta_formula"], {"H": -2, "O": -1})
        lost, = self.events("CH05-R004", "GRK", mass.calculate_mass(sequence="GK"), tags=["arg_coupling"])
        self.assertEqual(lost["product_structure"]["sequence"], ["G", "K"])
        observed = mass.calculate_mass(sequence="GK") + mass.calculate_mass(formula="C3H5NO")
        events = self.events("CH05-R015", "GK", observed, tags=["dcc_hosu"])
        self.assertEqual(len(events), 2)
        self.assertTrue(all("beta-Ala" in e["product_structure"]["sequence"] for e in events))

    def test_hoobt_cap_and_pyrrolidinamide(self):
        parent = mass.calculate_mass(sequence="GK")
        cap = mass.calculate_mass(formula="C7H3N3O")
        hit, = self.events("CH05-R017", "GK", parent + cap, tags=["dcc", "hoobt"])
        self.assertEqual(hit["site"], "N")
        self.assertFalse(self.events("CH05-R017", "GK", parent + cap, nterm="acetyl"))
        # C4H9N minus H2O is the shift relative to a free acid.
        pyr = mass.calculate_mass(formula="C4H9N") - mass.calculate_mass(formula="H2O")
        hit, = self.events("CH05-R032", "GK", parent + pyr,
                           tags=["phosphonium_coupling_reagent", "pyrrolidine_impurity"])
        self.assertEqual(hit["product_structure"]["cterm"], "pyrrolidinamide")
        self.assertFalse(self.events("CH05-R032", "GK", parent + pyr, cterm="amide"))

    def test_aspartimide_opening_isomers_and_precursor(self):
        parent = mass.calculate_mass(sequence="GDA")
        for rule, gain, reagent in [("CH06-R006", "C5H9N", "piperidine"),
                                     ("CH06-R007", "C5H9N", "piperidine"),
                                     ("CH06-R011", "CH2", "methanol"),
                                     ("CH06-R012", "CH2", "methanol")]:
            delta = mass.calculate_mass(formula=gain)
            if reagent == "piperidine":
                delta -= mass.calculate_mass(formula="O")
            event, = self.events(rule, "GDA", parent + delta, tags=[reagent, "base_exposure"])
            self.assertEqual(event["site"], [2, 3])
            self.assertFalse(self.events(rule, "GDP", mass.calculate_mass(sequence="GDP") + delta))
        for rule, sequence in [("CH06-R013", "GNA"), ("CH06-R018", "GQA")]:
            predicted = mass.calculate_mass(sequence=sequence) - mass.calculate_mass(formula="NH3")
            event, = self.events(rule, sequence, predicted, tags=["base_exposure"])
            self.assertEqual(event["site"], [2, 3])

    def test_n_terminal_cam_ring_distinguishes_acm(self):
        target = {"sequence": "CG", "nterm": "free", "cterm": "acid",
                  "modifications": [{"site": 1, "group": "Carbamidomethyl"}]}
        expected = mass.calculate_mass(sequence="CG") + mass.calculate_mass(formula="C2O")
        case = {"intended_product": target, "reaction": {"description": "bicarbonate",
                "stage_tags": ["ammonium_bicarbonate"]},
                "observations": [{"kind": "neutral_mono", "value": expected,
                                  "tolerance": 0.000003, "tolerance_unit": "Da"}]}
        output = diagnose(case, self.library)["observations"][0]
        event, = [e for h in output["matches"] for e in h["events"] if e["rule_id"] == "CH06-R041"]
        self.assertEqual(event["delta_formula"], {"H": -3, "N": -1})
        from diagnose_case import group_formula, formula_mass
        self.assertAlmostEqual(formula_mass(group_formula("Carbamidomethyl")), 57.02146372069, places=8)
        self.assertAlmostEqual(formula_mass(group_formula("Acm")), 71.03711378515, places=8)
        case["intended_product"] = dict(target, modifications=[{"site": 1, "group": "Acm"}])
        self.assertFalse(any(e["rule_id"] == "CH06-R041" for h in
            diagnose(case, self.library)["observations"][0]["matches"] for e in h["events"]))
