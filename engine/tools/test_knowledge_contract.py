"""Independent numbers, identity preservation, and derived-file integrity."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
from pyteomics import mass
import knowledge_source
from diagnose_case import load_rules, structure
from mass_catalog import group_formula

class KnowledgeContractTests(unittest.TestCase):
    def test_component_identities_are_not_folded_by_mass(self):
        self.assertEqual(group_formula('Carbamidomethyl'),{'C':2,'H':3,'N':1,'O':1})
        self.assertEqual(group_formula('Acm'),{'C':3,'H':5,'N':1,'O':1})
        self.assertEqual(group_formula('Acm'),knowledge_source.components()['residues']['A']['residue_formula'])
        self.assertNotEqual(knowledge_source.components()['covalent_groups']['Acm']['component_kind'],
                            knowledge_source.components()['residues']['A']['component_kind'])

    def test_formula_mass_and_termini_are_independent(self):
        from mass_ledger import formula_mass
        f=structure({'sequence':'C','nterm':'free','cterm':'amide','modifications':[{'site':1,'group':'Carbamidomethyl'}]})[1]
        # Independent isotope tables differ slightly for sulfur; sub-microdalton tolerance.
        self.assertAlmostEqual(formula_mass(f),mass.calculate_mass(sequence='C')+57.02146372069-0.98401558291,places=6)

    def test_frozen_runtime_preserves_every_legacy_and_branch_id(self):
        bundle,_=knowledge_source.load_runtime()
        self.assertEqual(len(bundle['rules']),450)
        self.assertEqual(sum(bool(r['template']) for r in bundle['rules']),240)
        self.assertEqual(len({b['branch_id'] for r in bundle['rules'] for b in r['branches']}),308)
        ids={r['rule_id'] for r in bundle['rules']}
        self.assertTrue({'CH09-R024','SUP-056','APP-I-R004'}<=ids)
        rules,_=load_rules();self.assertEqual(rules['SUP-056']['alias_of'],'APP-I-R004')
        facts={r['rule_id']:r for r in bundle['appendix_dispositions']}
        self.assertEqual(facts['AI-040']['status'],'compiled_fact')
        self.assertEqual(facts['AI-081']['record_kind'],'component_depiction')

    def test_tampered_export_and_stale_source_are_detected(self):
        root=knowledge_source.ROOT
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)
            shutil.copytree(root/'knowledge',target/'knowledge')
            shutil.copytree(root/'runtime',target/'runtime')
            manifest=json.loads((target/'runtime/manifest.json').read_text(encoding='utf-8'))
            for name in manifest['compatibility_exports']:
                dst=target/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/name,dst)
            with patch.object(knowledge_source,'ROOT',target):
                knowledge_source.load_runtime()
                p=target/'SUPPLEMENTS.json';original=p.read_bytes();p.write_bytes(original+b' ')
                with self.assertRaisesRegex(ValueError,'Modified derived export'):knowledge_source.load_runtime()
                p.write_bytes(original)
                q=target/'knowledge/components.json';q.write_bytes(q.read_bytes()+b' ')
                with self.assertRaisesRegex(ValueError,'Stale knowledge runtime'):knowledge_source.load_runtime()
