import unittest
from pyteomics import mass
from diagnose_case import diagnose, load_rules, structure


class ChapterThirteenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = load_rules()

    @staticmethod
    def target(seq="GCCG", kind="disulfide", mods=()):
        return dict(sequence=seq, nterm="free", cterm="acid", modifications=list(mods),
                    **({"bonds":[{"kind":kind,"sites":[2,3]}],"topology":"crosslinked"} if kind else {}))

    def result(self, target, value, tolerance=0.000004):
        return diagnose({"intended_product":target,"reaction":{"description":"CH13 explicit sulfur topology"},
                         "observations":[{"kind":"neutral_mono","value":value,"tolerance":tolerance,"tolerance_unit":"Da"}]},self.library)

    @staticmethod
    def events(result,n):
        return [e for h in result["observations"][0]["matches"] for e in h["events"] if e["rule_id"]==f"CH13-R{n:03}"]

    def test_four_cys_networks_match_independent_three_pairing_answers(self):
        target = self.target("CCCC")
        target["bonds"] = [{"kind":"disulfide","sites":[1,2]},{"kind":"disulfide","sites":[3,4]}]
        observed = mass.calculate_mass(sequence="CCCC")-2*mass.calculate_mass(formula="H2")
        found = self.events(self.result(target,observed),1)
        actual = {frozenset(tuple(b["sites"]) for b in e["product_structure"]["bonds"]) for e in found}
        self.assertEqual(actual,{frozenset({(1,3),(2,4)}),frozenset({(1,4),(2,3)})})
        self.assertTrue(all(e["delta_formula"]=={} and e["elementary_exchange_steps"]=="not_enumerated" for e in found))
        blocked = self.target("CCCG",mods=[{"site":3,"group":"Acm"}])
        blocked["bonds"][0]["sites"]=[1,2]
        self.assertFalse(self.events(self.result(blocked,mass.calculate_mass(sequence="CCCG")+mass.calculate_mass(formula="C3H5NO")-mass.calculate_mass(formula="H2")),1))

    def test_network_budget_rejects_without_silent_truncation(self):
        target = self.target("C"*16)
        target["bonds"]=[{"kind":"disulfide","sites":[i,i+1]} for i in range(1,17,2)]
        with self.assertRaisesRegex(ValueError,"2027025 states"):
            self.result(target,mass.calculate_mass(sequence="C"*16)-8*mass.calculate_mass(formula="H2"))

    def test_disulfide_hydrolysis_beta_and_alpha_elimination_use_whole_backbone(self):
        free = mass.calculate_mass(sequence="GCCG")
        target = self.target()
        for n,observed,delta in [(4,free+mass.calculate_mass(formula="O"),{"H":2,"O":1}),
                                  (6,free-mass.calculate_mass(formula="H2"),{}),
                                  (9,free-mass.calculate_mass(formula="H2"),{})]:
            with self.subTest(rule=n):
                found = self.events(self.result(target,observed),n)
                self.assertEqual(len(found),2)
                self.assertTrue(all(e["delta_formula"]==delta for e in found))
                self.assertTrue(all(e["product_structure"]["sequence"] or n==6 for e in found))
        target["alpha_stereochemistry"]={"2":"L","3":"D"}
        found = self.events(self.result(target,free-mass.calculate_mass(formula="H2")),6)
        for e in found:
            structure(e["product_structure"])
            dha=e["product_structure"]["sequence"].index("Dha")+1
            self.assertNotIn(str(dha),e["product_structure"].get("alpha_stereochemistry",{}))

    def test_named_persulfide_and_thioaldehyde_fates_preserve_actual_references(self):
        parent=mass.calculate_mass(sequence="GCG")
        for n,group,formula,delta in [(7,"Persulfide","",{"S":-1}),
                                      (8,"Persulfide","O",{"S":-1,"O":1}),
                                      (10,"CysThioaldehyde","H-2OS-1",{"S":-1,"O":1})]:
            with self.subTest(rule=n):
                observed=parent+(mass.calculate_mass(formula=formula) if formula else 0)
                event,=self.events(self.result(self.target("GCG",None,[{"site":2,"group":group}]),observed),n)
                self.assertEqual(event["delta_formula"],delta)
                self.assertFalse(self.events(self.result(self.target("GCG",None),observed),n))
        target=self.target(kind=None,mods=[{"site":2,"group":"CysSulfenic"},{"site":3,"group":"CysSulfenic"}])
        found=self.events(self.result(target,mass.calculate_mass(sequence="GCCG")+mass.calculate_mass(formula="O2")),5)
        self.assertEqual(len(found),2)
        self.assertTrue(all(e["site"]==[2,3] and e["delta_formula"]=={} for e in found))

    def test_inserted_sulfur_and_reduction_state_mass_balance(self):
        free=mass.calculate_mass(sequence="GCCG")
        ss=free-mass.calculate_mass(formula="H2")
        for n,delta in [(11,"S"),(13,"S2"),(14,"S")]:
            matched=self.events(self.result(self.target(),ss+mass.calculate_mass(formula=delta)),n)
            self.assertEqual(len(matched),2 if n==14 else 1)
            event=matched[0]
            self.assertEqual(event["delta_formula"],{"S":2 if n==13 else 1})
            structure(event["product_structure"])
        event,=self.events(self.result(self.target(kind="trisulfide"),free),15)
        self.assertEqual(event["delta_formula"],{"H":2,"S":-1})
        self.assertNotIn("bonds",event["product_structure"])
        lanth=self.events(self.result(self.target(),ss-mass.calculate_mass(formula="S")),16)
        self.assertEqual(len(lanth),2)
        self.assertTrue(all(e["delta_formula"]=={"S":-1} for e in lanth))

    def test_low_resolution_retains_sulfur_and_two_oxygen_isobars(self):
        ss=mass.calculate_mass(sequence="GCCG")-mass.calculate_mass(formula="H2")
        observed=ss+mass.calculate_mass(formula="S")
        broad=self.result(self.target(),observed,tolerance=0.02)
        ids={e["rule_id"] for h in broad["observations"][0]["matches"] for e in h["events"]}
        self.assertIn("CH13-R011",ids)
        self.assertIn("CH09-R007",ids)
        narrow=self.result(self.target(),observed)
        ids={e["rule_id"] for h in narrow["observations"][0]["matches"] for e in h["events"]}
        self.assertIn("CH13-R011",ids)
        self.assertNotIn("CH09-R007",ids)

    def test_dha_amine_michael_addition_records_unsaturation_and_chain_count(self):
        for n,seq,formula in [(17,["G","Dha","Orn","G"],"C12H21N5O5"),
                               (18,["G","Dha","K","G"],"C13H23N5O5")]:
            parent=mass.calculate_mass(formula=formula)
            for observed,chains in [(parent,1),(2*parent,2)]:
                with self.subTest(rule=n,chains=chains):
                    event,=self.events(self.result(self.target(seq,None),observed),n)
                    self.assertEqual(len(event["product_structure"]["chains"]),chains)
                    self.assertEqual(event["product_structure"]["bonds"][0]["consumed_unsaturation"]["site"],2)
                    self.assertEqual(event["product_structure"]["bonds"][0]["new_alpha_configuration"],"not_assigned")


if __name__=="__main__":
    unittest.main()
