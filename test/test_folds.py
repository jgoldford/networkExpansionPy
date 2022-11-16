import unittest
import networkExpansionPy.lib as ne
import networkExpansionPy.folds as nf
import pandas as pd
from scipy.sparse import csr_matrix
from pprint import pprint
# from pandas.testing import assert_frame_equal

class TestIndependentFunctions(unittest.TestCase):

    maxDiff = None

    def test_rule2rn(self):
        rn2rules = {
            'R0': {frozenset({'F0'}), frozenset({'F10'}), frozenset({'F11'}), frozenset({'F12'})},
            'R1': {frozenset({'F1'})},
            'R2': {frozenset({'F11'}), frozenset({'F12'}), frozenset({'F2'})},
            'R3': {frozenset({'F3'})},
            'R4': {frozenset({'F4'})},
            'R5': {frozenset({'F5'})},
            'R6': {frozenset({'F6'})},
            'R7': {frozenset({'F7'})},
            'R8': {frozenset({'F8'})},
            'R9': {frozenset({'F9'})}
            }

        expected_rule2rns = {
            frozenset({'F12'}): {'R2', 'R0'},
            frozenset({'F11'}): {'R2', 'R0'},
            frozenset({'F10'}): {'R0'},
            frozenset({'F0'}): {'R0'},
            frozenset({'F1'}): {'R1'},
            frozenset({'F2'}): {'R2'},
            frozenset({'F3'}): {'R3'},
            frozenset({'F4'}): {'R4'},
            frozenset({'F5'}): {'R5'},
            frozenset({'F6'}): {'R6'},
            frozenset({'F7'}): {'R7'},
            frozenset({'F8'}): {'R8'},
            frozenset({'F9'}): {'R9'}
            }

        self.assertEqual(expected_rule2rns, nf.rule2rn(rn2rules))

    def test_subset_rule2rn(self):
        
        rule2rn = {
            frozenset({'F12'}): {'R2', 'R0'},
            frozenset({'F11'}): {'R2', 'R0'},
            frozenset({'F10','F11'}): {'R0'},
            frozenset({'F0'}): {'R0'},
            frozenset({'F1'}): {'R1'},
            frozenset({'F2'}): {'R2'},
            frozenset({'F3'}): {'R3'},
            frozenset({'F4'}): {'R4'},
            frozenset({'F5'}): {'R5'},
            frozenset({'F6'}): {'R6'},
            frozenset({'F7'}): {'R7'},
            frozenset({'F8'}): {'R8'},
            frozenset({'F9'}): {'R9'}
            }


        folds = {'F10'}
        expected = {}
        self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

        folds = {'F10','F11'}
        expected = {
            frozenset({'F11'}): {'R2', 'R0'}, 
            frozenset({'F10','F11'}): {'R0'}
            }
        self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

        folds = {'F9','F11'}
        expected = {
            frozenset({'F11'}): {'R2', 'R0'}, 
            frozenset({'F9'}): {'R9'}
            }
        self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

    def test_rule2nextrns(self):

        rule2rn = {
            frozenset({'F12'}): {'R2', 'R0'},
            frozenset({'F11'}): {'R2', 'R0'},
            frozenset({'F10','F11'}): {'R0'},
            frozenset({'F0'}): {'R0'},
            frozenset({'F1'}): {'R1'},
            frozenset({'F2'}): {'R2'},
            frozenset({'F3'}): {'R3'},
            frozenset({'F4'}): {'R4'},
            frozenset({'F5'}): {'R5'},
            frozenset({'F6'}): {'R6'},
            frozenset({'F7'}): {'R7'},
            frozenset({'F8'}): {'R8'},
            frozenset({'F9'}): {'R9'}
            }

        folds = {'F0'}
        expected = {
            frozenset({'F12'}): {'R2', 'R0'},
            frozenset({'F11'}): {'R2', 'R0'},
            frozenset({'F1'}): {'R1','R0'},
            frozenset({'F2'}): {'R2','R0'},
            frozenset({'F3'}): {'R3','R0'},
            frozenset({'F4'}): {'R4','R0'},
            frozenset({'F5'}): {'R5','R0'},
            frozenset({'F6'}): {'R6','R0'},
            frozenset({'F7'}): {'R7','R0'},
            frozenset({'F8'}): {'R8','R0'},
            frozenset({'F9'}): {'R9','R0'}
            }

        self.assertEqual(expected, nf.rule2nextrns(folds, rule2rn))

        folds = {'F0', 'F1','F2','F3','F4','F5'}
        expected = {
            frozenset({'F6'}): {'R6','R0','R1','R2','R3','R4','R5'},
            frozenset({'F7'}): {'R7','R0','R1','R2','R3','R4','R5'},
            frozenset({'F8'}): {'R8','R0','R1','R2','R3','R4','R5'},
            frozenset({'F9'}): {'R9','R0','R1','R2','R3','R4','R5'}
            }

        self.assertEqual(expected, nf.rule2nextrns(folds, rule2rn))

    def test_create_equal_rule_groups(self):
        pass


# class TestGlobalFoldNetworkIrreversible(unittest.TestCase):

#     maxDiff = None ## allows full output of failed test differences

#     def setUp(self):
#         reactions = 10
#         compounds = 11
#         rids = ['R' + str(x) for x in range(reactions)]
#         cids = ['C' + str(x) for x in range(compounds)]
#         folds = ['F' + str(x) for x in range(reactions)]
#         network = {'rn':[],'direction':[],'cid':[],'s':[]}
#         fold_rules  = {'rn': [],'rule':[]}
#         i = 0
#         for r in rids:
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i])    
#             network['s'].append(-1)
            
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i+1])    
#             network['s'].append(1)
#             fold_rules['rn'].append(r)
#             fold_rules['rule'].append(folds[i])
            
#             i = i +1

#         self.network = pd.DataFrame(network)
#         self.fold_rules = pd.DataFrame(fold_rules)
#         self.fold_rules['fold_sets'] = self.fold_rules.rule.apply(lambda x: set(x.split('_')))
#         self.rn2rules = {d["rn"]:{frozenset(d["fold_sets"])} for d in self.fold_rules.to_dict(orient="records")}

#         ## Create Metabolism
#         self.met = ne.GlobalMetabolicNetwork(metabolism="dev")
#         self.met.network = self.network

#     def test_GlobalFoldNetwork_create_foldrules2rn(self):
#         fold_independent_rns = set()
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         expected_rule2rns = {frozenset({'F0'}): {'R0'},
#                         frozenset({'F1'}): {'R1'},
#                         frozenset({'F2'}): {'R2'},
#                         frozenset({'F3'}): {'R3'},
#                         frozenset({'F4'}): {'R4'},
#                         frozenset({'F5'}): {'R5'},
#                         frozenset({'F6'}): {'R6'},
#                         frozenset({'F7'}): {'R7'},
#                         frozenset({'F8'}): {'R8'},
#                         frozenset({'F9'}): {'R9'}}
#         self.assertEqual(foldnet.rule2rns, expected_rule2rns)

#     def test_FoldMetabolism_rule_order_C0_no_indepdendent(self):
#         fold_independent_rns = set()
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
#                                 frozenset({'F1'}): {'R1'},
#                                 frozenset({'F2'}): {'R2'},
#                                 frozenset({'F3'}): {'R3'},
#                                 frozenset({'F4'}): {'R4'},
#                                 frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
#                                 'R1': {frozenset({'F1'})},
#                                 'R2': {frozenset({'F2'})},
#                                 'R3': {frozenset({'F3'})},
#                                 'R4': {frozenset({'F4'})},
#                                 'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C0': 0,
#                                 'C1': 2,
#                                 'C2': 3,
#                                 'C3': 4,
#                                 'C4': 5,
#                                 'C5': 6,
#                                 'C6': 7,
#                                 'C7': 8,
#                                 'C8': 9,
#                                 'C9': 10,
#                                 'C10': 11},
#                                 'rns': {'R0': 2,
#                                 'R1': 3,
#                                 'R2': 4,
#                                 'R3': 5,
#                                 'R4': 6,
#                                 'R5': 7,
#                                 'R6': 8,
#                                 'R7': 9,
#                                 'R8': 10,
#                                 'R9': 11},
#                                 'folds': {'fold_independent': 0,
#                                 'F0': 2,
#                                 'F1': 3,
#                                 'F2': 4,
#                                 'F3': 5,
#                                 'F4': 6,
#                                 'F5': 7,
#                                 'F6': 8,
#                                 'F7': 9,
#                                 'F8': 10,
#                                 'F9': 11}}

#         self.assertEqual(iteration_dict, expected_iteration_dict)
#         self.assertEqual(current, expected_current)
        
#     def test_FoldMetabolism_rule_order_C0_independent_R0R1(self):
#         fold_independent_rns = set(["R0","R1"])
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
#                                 frozenset({'F1'}): {'R1'},
#                                 frozenset({'F2'}): {'R2'},
#                                 frozenset({'F3'}): {'R3'},
#                                 frozenset({'F4'}): {'R4'},
#                                 frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
#                                 'R1': {frozenset({'F1'})},
#                                 'R2': {frozenset({'F2'})},
#                                 'R3': {frozenset({'F3'})},
#                                 'R4': {frozenset({'F4'})},
#                                 'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C0': 0,
#                                 'C1': 1,
#                                 'C2': 1,
#                                 'C3': 2,
#                                 'C4': 3,
#                                 'C5': 4,
#                                 'C6': 5,
#                                 'C7': 6,
#                                 'C8': 7,
#                                 'C9': 8,
#                                 'C10': 9},
#                                 'rns': {'R0': 1,
#                                 'R1': 1,
#                                 'R2': 2,
#                                 'R3': 3,
#                                 'R4': 4,
#                                 'R5': 5,
#                                 'R6': 6,
#                                 'R7': 7,
#                                 'R8': 8,
#                                 'R9': 9},
#                                 'folds': {'fold_independent': 0,
#                                 'F0': 2,
#                                 'F1': 3,
#                                 'F2': 4,
#                                 'F3': 5,
#                                 'F4': 6,
#                                 'F5': 7,
#                                 'F6': 8,
#                                 'F7': 9,
#                                 'F8': 10,
#                                 'F9': 11}}

#     def test_FoldMetabolism_rule_order_C0_independent_R3R5(self):
#         fold_independent_rns = set(["R3","R5"])
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
#                                 frozenset({'F1'}): {'R1'},
#                                 frozenset({'F2'}): {'R2'},
#                                 frozenset({'F3'}): {'R3'},
#                                 frozenset({'F4'}): {'R4'},
#                                 frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
#                                 'R1': {frozenset({'F1'})},
#                                 'R2': {frozenset({'F2'})},
#                                 'R3': {frozenset({'F3'})},
#                                 'R4': {frozenset({'F4'})},
#                                 'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C0': 0,
#                                 'C1': 2,
#                                 'C2': 3,
#                                 'C3': 4,
#                                 'C4': 4,
#                                 'C5': 5,
#                                 'C6': 5,
#                                 'C7': 6,
#                                 'C8': 7,
#                                 'C9': 8,
#                                 'C10': 9},
#                                 'rns': {'R0': 2,
#                                 'R1': 3,
#                                 'R2': 4,
#                                 'R3': 4,
#                                 'R4': 5,
#                                 'R5': 5,
#                                 'R6': 6,
#                                 'R7': 7,
#                                 'R8': 8,
#                                 'R9': 9},
#                                 'folds': {'fold_independent': 0,
#                                 'F0': 2,
#                                 'F1': 3,
#                                 'F2': 4,
#                                 'F3': 5,
#                                 'F4': 6,
#                                 'F5': 7,
#                                 'F6': 8,
#                                 'F7': 9,
#                                 'F8': 10,
#                                 'F9': 11}}

#     def test_FoldMetabolism_rule_order_C5_no_independent(self):
#         fold_independent_rns = set()
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C5'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C10', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C5': 0, 'C6': 2, 'C7': 3, 'C8': 4, 'C9': 5, 'C10': 6},
#                             'rns': {'R5': 2, 'R6': 3, 'R7': 4, 'R8': 5, 'R9': 6},
#                             'folds': {'fold_independent': 0, 'F5': 2, 'F6': 3, 'F7': 4, 'F8': 5, 'F9': 6}}
#         self.assertEqual(iteration_dict, expected_iteration_dict)
#         self.assertEqual(current, expected_current)

# class TestGlobalFoldNetworkReversible(unittest.TestCase):

#     def setUp(self):
#         reactions = 10
#         compounds = 11
#         rids = ['R' + str(x) for x in range(reactions)]
#         cids = ['C' + str(x) for x in range(compounds)]
#         folds = ['F' + str(x) for x in range(reactions)]
#         network = {'rn':[],'direction':[],'cid':[],'s':[]}
#         fold_rules  = {'rn': [],'rule':[]}
#         i = 0
#         for r in rids:
#             ## Forward i
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i])    
#             network['s'].append(-1)
#             ## Backwards i
#             network['rn'].append(r)
#             network['direction'].append('backward')
#             network['cid'].append(cids[i])    
#             network['s'].append(1)
            
#             ## Forward i+1
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i+1])    
#             network['s'].append(1)
#             fold_rules['rn'].append(r)
#             fold_rules['rule'].append(folds[i])
#             ## Backwards i+1
#             network['rn'].append(r)
#             network['direction'].append('backward')
#             network['cid'].append(cids[i+1])    
#             network['s'].append(-1)
#             fold_rules['rn'].append(r)
#             fold_rules['rule'].append(folds[i])
            
#             i = i +1

#         self.network = pd.DataFrame(network)
#         self.fold_rules = pd.DataFrame(fold_rules)
#         self.fold_rules['fold_sets'] = self.fold_rules.rule.apply(lambda x: set(x.split('_')))
#         self.rn2rules = {d["rn"]:{frozenset(d["fold_sets"])} for d in self.fold_rules.to_dict(orient="records")}

#         ## Create Metabolism
#         self.met = ne.GlobalMetabolicNetwork(metabolism="dev")
#         self.met.network = self.network

#     def test_FoldMetabolism_rule_order_C0_reversible(self):
#         fold_independent_rns = set()
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
#                                 frozenset({'F1'}): {'R1'},
#                                 frozenset({'F2'}): {'R2'},
#                                 frozenset({'F3'}): {'R3'},
#                                 frozenset({'F4'}): {'R4'},
#                                 frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
#                                 'R1': {frozenset({'F1'})},
#                                 'R2': {frozenset({'F2'})},
#                                 'R3': {frozenset({'F3'})},
#                                 'R4': {frozenset({'F4'})},
#                                 'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C0': 0,
#                                 'C1': 2,
#                                 'C2': 3,
#                                 'C3': 4,
#                                 'C4': 5,
#                                 'C5': 6,
#                                 'C6': 7,
#                                 'C7': 8,
#                                 'C8': 9,
#                                 'C9': 10,
#                                 'C10': 11},
#                                 'rns': {'R0': 2,
#                                 'R1': 3,
#                                 'R2': 4,
#                                 'R3': 5,
#                                 'R4': 6,
#                                 'R5': 7,
#                                 'R6': 8,
#                                 'R7': 9,
#                                 'R8': 10,
#                                 'R9': 11},
#                                 'folds': {'fold_independent': 0,
#                                 'F0': 2,
#                                 'F1': 3,
#                                 'F2': 4,
#                                 'F3': 5,
#                                 'F4': 6,
#                                 'F5': 7,
#                                 'F6': 8,
#                                 'F7': 9,
#                                 'F8': 10,
#                                 'F9': 11}}

#         self.assertEqual(iteration_dict, expected_iteration_dict)
#         self.assertEqual(current, expected_current)

#     def test_FoldMetabolism_rule_order_C5_reversible(self):
#         fold_independent_rns = set([])
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C5'])
#         fm.seed_folds = set([])

#         expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
#                                 frozenset({'F1'}): {'R1'},
#                                 frozenset({'F2'}): {'R2'},
#                                 frozenset({'F3'}): {'R3'},
#                                 frozenset({'F4'}): {'R4'},
#                                 frozenset({'F5'}): {'R5'},
#                                 frozenset({'F6'}): {'R6'},
#                                 frozenset({'F7'}): {'R7'},
#                                 frozenset({'F8'}): {'R8'},
#                                 frozenset({'F9'}): {'R9'}}
#         self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

#         expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
#                                 'R1': {frozenset({'F1'})},
#                                 'R2': {frozenset({'F2'})},
#                                 'R3': {frozenset({'F3'})},
#                                 'R4': {frozenset({'F4'})},
#                                 'R5': {frozenset({'F5'})},
#                                 'R6': {frozenset({'F6'})},
#                                 'R7': {frozenset({'F7'})},
#                                 'R8': {frozenset({'F8'})},
#                                 'R9': {frozenset({'F9'})}}
#         self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

#         current, iteration_dict, metadict = fm.rule_order()
#         expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
#                             'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
#                             'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

#         expected_iteration_dict = {'cpds': {'C0': 6,
#                                 'C1': 5,
#                                 'C2': 4,
#                                 'C3': 3,
#                                 'C4': 2,
#                                 'C5': 0,
#                                 'C6': 7,
#                                 'C7': 8,
#                                 'C8': 9,
#                                 'C9': 10,
#                                 'C10': 11},
#                                 'rns': {'R0': 6,
#                                 'R1': 5,
#                                 'R2': 4,
#                                 'R3': 3,
#                                 'R4': 2,
#                                 'R5': 7,
#                                 'R6': 8,
#                                 'R7': 9,
#                                 'R8': 10,
#                                 'R9': 11},
#                                 'folds': {'fold_independent': 0,
#                                 'F0': 2,
#                                 'F1': 3,
#                                 'F2': 4,
#                                 'F3': 5,
#                                 'F4': 6,
#                                 'F5': 7,
#                                 'F6': 8,
#                                 'F7': 9,
#                                 'F8': 10,
#                                 'F9': 11}}

# class TestGlobalFoldNetworkTwoFoldsSimultaneouslyNeeded(unittest.TestCase):

#     def setUp(self):
#         reactions = 10
#         compounds = 11
#         rids = ['R' + str(x) for x in range(reactions)]
#         cids = ['C' + str(x) for x in range(compounds)]
#         folds = ['F' + str(x) for x in range(reactions)]
#         network = {'rn':[],'direction':[],'cid':[],'s':[]}
#         fold_rules  = {'rn': [],'rule':[]}
#         i = 0
#         for r in rids:
#             ## Forward i
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i])    
#             network['s'].append(-1)

#             ## Forward i+1
#             network['rn'].append(r)
#             network['direction'].append('forward')
#             network['cid'].append(cids[i+1])    
#             network['s'].append(1)
#             fold_rules['rn'].append(r)
#             fold_rules['rule'].append(folds[i])
            
#             i = i +1

#         self.network = pd.DataFrame(network)
#         self.fold_rules = pd.DataFrame(fold_rules)
#         self.fold_rules['fold_sets'] = self.fold_rules.rule.apply(lambda x: set(x.split('_')))
#         self.rn2rules = {d["rn"]:{frozenset(d["fold_sets"])} for d in self.fold_rules.to_dict(orient="records")}
#         self.rn2rules["R0"] = {frozenset({'F0','F10'})}

#         ## Create Metabolism
#         self.met = ne.GlobalMetabolicNetwork(metabolism="dev")
#         self.met.network = self.network

#     def test_FoldMetabolism_rule_order_R0_needs_2_folds(self):
#         fold_independent_rns = set([])
#         foldnet = nf.GlobalFoldNetwork(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldnet)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#     def test_exit_at_correct_iteration_1(self):
#         pass