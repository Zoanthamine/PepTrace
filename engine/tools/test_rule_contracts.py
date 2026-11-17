import unittest
from copy import deepcopy
from knowledge_source import load_runtime
from rule_contracts import validate_contract

class RuleContractTests(unittest.TestCase):
    def test_every_compiled_branch_and_legacy_id_has_a_route(self):
        bundle,_=load_runtime();contracts=[b['contract'] for r in bundle['rules'] for b in r['branches']]
        self.assertEqual(len(contracts),308)
        for c in contracts:validate_contract(c)
        self.assertEqual(len([r for r in bundle['rules'] if r['execution_route']]),450)
    def test_unknown_condition_and_empty_consumption_are_rejected(self):
        contract=next(b['contract'] for r in load_runtime()[0]['rules'] for b in r['branches'])
        for edit in [{'conditions':{**contract['conditions'],'unknown_policy':'drop'}},{'consumed_resources':[]}]:
            with self.subTest(edit=edit),self.assertRaises(ValueError):validate_contract({**deepcopy(contract),**edit})
