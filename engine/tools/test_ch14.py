import unittest
from pyteomics import mass
from diagnose_case import diagnose,load_rules,structure
from carbonyl_components import normalize_carbonyl_donors


class ChapterFourteenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library=load_rules()

    @staticmethod
    def target(seq,mods=(),cterm="acid"):
        return dict(sequence=seq,nterm="free",cterm=cterm,modifications=list(mods))

    def result(self,target,value,tags=(),donors=None,start=None,observation=None):
        reaction=dict(description="CH14 independent model",stage_tags=list(tags))
        if donors is not None:
            reaction["carbonyl_donors"]=donors
        case={"intended_product":target,"reaction":reaction,"observations":[observation or
              {"kind":"neutral_mono","value":value,"tolerance":0.000004,"tolerance_unit":"Da"}]}
        if start is not None:
            case["starting_material"]=start
        return diagnose(case,self.library)

    @staticmethod
    def events(result,n,bucket="matches"):
        return [e for h in result["observations"][0][bucket] for e in h["events"] if e["rule_id"]==f"CH14-R{n:03}"]

    def test_fmoc_retention_derives_explicit_starting_protection(self):
        target=self.target("GK")
        observed=mass.calculate_mass(sequence="GK")+mass.calculate_mass(formula="C15H10O2")
        start=self.target("GK",[{"site":"N","group":"Fmoc"}])
        event,=self.events(self.result(target,observed,["dcm","piperidine","fmoc_deprotection"],start=start),2)
        self.assertEqual(event["condition"],"compatible")
        self.assertEqual(event["product_structure"]["modifications"][0]["group"],"Fmoc")
        result=self.result(target,observed,["dcm","piperidine","fmoc_deprotection"],start=target)
        self.assertFalse(self.events(result,2))
        event,=self.events(result,2,"condition_conflicts")
        self.assertIn("n_fmoc_precursor:N",event["conflicting_conditions"])

    def test_nca_chloride_branch_has_combined_terminal_state_and_independent_formula(self):
        for residue in "GAP":
            with self.subTest(residue=residue):
                parent=mass.calculate_mass(sequence=residue)
                target=self.target(residue,[{"site":"N","group":"NCA"}])
                observed=parent+mass.calculate_mass(formula="CH-2Cl2")
                event,=self.events(self.result(target,observed),6)
                self.assertEqual(event["delta_formula"],{"Cl":2,"O":-1})
                self.assertEqual(event["product_structure"]["modifications"][0]["group"],"NChloroformylAcidChloride")
                self.assertAlmostEqual(mass.calculate_mass(composition=structure(target)[1]),
                                       parent+mass.calculate_mass(formula="CH-2O"),places=7)
                structure(event["product_structure"])
        with self.assertRaisesRegex(ValueError,"one alpha amino acid"):
            structure(self.target("GG",[{"site":"N","group":"NCA"}]))
        with self.assertRaisesRegex(ValueError,"no separate C-terminal state"):
            structure(self.target("G",[{"site":"N","group":"NCA"},{"site":"C","group":"C-methyl-ester"}]))

    def test_ester_and_reductive_alkylation_sites_keep_amine_hydrogens(self):
        seq="GDEK"
        parent=mass.calculate_mass(sequence=seq)
        found=self.events(self.result(self.target(seq),parent+mass.calculate_mass(formula="CH2")),7)
        self.assertEqual({e["site"] for e in found},{2,3,"C"})
        for n,formula in [(8,"CH2"),(9,"C2H4"),(10,"C2H4"),(11,"C7H6")]:
            found=self.events(self.result(self.target(seq),parent+mass.calculate_mass(formula=formula)),n)
            self.assertEqual({e["site"] for e in found},{"N",4})
        self.assertFalse(self.events(self.result(self.target("PG"),mass.calculate_mass(sequence="PG")+mass.calculate_mass(formula="C2H4")),9))

    def test_acetone_his_imine_and_cyclic_amide_are_distinct(self):
        parent=mass.calculate_mass(sequence="HG")
        event,=self.events(self.result(self.target("HG"),parent+mass.calculate_mass(formula="C3H4")),15)
        self.assertEqual(event["site"],1)
        events=self.events(self.result(self.target("HG"),parent+mass.calculate_mass(formula="C3H4")),16)
        self.assertEqual(len(events),2)
        self.assertEqual({e["product_structure"]["topology"] for e in events},{"linear","local_cycle"})
        amide=self.events(self.result(self.target("G",cterm="amide"),mass.calculate_mass(formula="C5H10N2O")),16)
        self.assertEqual(len(amide),2)
        self.assertTrue(any(e["site"]==["N",1,"C"] for e in amide))
        gp=self.events(self.result(self.target("GP"),mass.calculate_mass(sequence="GP")+mass.calculate_mass(formula="C3H4")),16)
        self.assertTrue(all(e["product_structure"]["topology"]=="linear" for e in gp))

    def test_aminooxy_carbonyl_condensation_validates_structure_and_ez(self):
        target=self.target(["Aoa","G"])
        parent=mass.calculate_mass(formula="C4H8N2O4")
        event,=self.events(self.result(target,parent+mass.calculate_mass(formula="C3H4")),17)
        self.assertEqual(event["site"],"N")
        for smiles,formula,count in [("CC(C)=O","C3H4",1),("CC=O","C2H2",2),("C=O","C",1)]:
            with self.subTest(smiles=smiles):
                observed=parent+mass.calculate_mass(formula=formula)
                result=self.result(target,observed,donors=[{"smiles":smiles}])
                found=self.events(result,18)
                self.assertEqual(len(found),count)
                if count==2:
                    self.assertEqual({e["product_structure"]["condensation"]["oxime_configuration"] for e in found},{"E","Z"})
        missing=self.result(target,parent+mass.calculate_mass(formula="C2H2"))
        self.assertIn("CH14-R018",missing["scope"]["carbonyl_component_rules_skipped"])
        for smiles in ["CC(=O)O","NC=O","CC(C)=O.C=O","[13CH3]C(C)=O"]:
            with self.subTest(invalid=smiles),self.assertRaises(ValueError):
                normalize_carbonyl_donors([{"smiles":smiles}])
        self.assertEqual(normalize_carbonyl_donors([{"smiles":"CC=O"}]),
                         normalize_carbonyl_donors([{"smiles":"[H]C(=O)C"}]))

    def test_mtbe_met_sulfonium_uses_intrinsic_charge_instead_of_neutral_tbu_delta(self):
        parent=mass.calculate_mass(sequence="GWYCMG")
        events=self.events(self.result(self.target("GWYCMG"),parent+mass.calculate_mass(formula="C4H8")),20)
        self.assertEqual({e["site"] for e in events},{2,3,4})
        parent=mass.calculate_mass(sequence="GMG")
        observation={"kind":"mz","value":parent+mass.calculate_mass(formula="C4H9")-0.000548579909,
                     "charge":1,"intrinsic_charge":1,"ions":{},"tolerance":0.000004,"tolerance_unit":"Da"}
        event,=self.events(self.result(self.target("GMG"),None,observation=observation),20)
        self.assertEqual(event["delta_formula"],{"C":4,"H":9})
        self.assertEqual(event["intrinsic_charge"],1)
        self.assertFalse(self.events(self.result(self.target("GMG"),parent+mass.calculate_mass(formula="C4H8")),20))

    def test_tfe_acylation_and_unconfirmed_nominal_leads_remain_separate(self):
        parent=mass.calculate_mass(sequence="GK")
        observed=parent+mass.calculate_mass(formula="C2H-1F3O")
        event,=self.events(self.result(self.target("GK"),observed),21)
        self.assertEqual(event["site"],2)
        for value,rule in [(27,"CH14-R012"),(26,"CH14-R022")]:
            result=diagnose({"reaction":{"description":"book unknown observation"},
                             "observations":[{"kind":"nominal_delta","value":value,"reference":"intended_product"}]},self.library)
            leads=[l for l in result["observations"][0]["reference_only_delta_leads"] if l["rule_id"]==rule]
            self.assertEqual(len(leads),1)
            self.assertIsNone(leads[0]["theoretical_delta_mono"])


if __name__=="__main__":
    unittest.main()
