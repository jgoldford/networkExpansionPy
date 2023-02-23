import unittest
import networkExpansionPy.lib as ne
import networkExpansionPy.folds as nf
import pandas as pd
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
from pprint import pprint
import random

import warnings
# warnings.filterwarnings('ignore', category=SparseEfficiencyWarning) ## must be applied within each test
unittest.defaultTestLoader.sortTestMethodsUsing = lambda *args: -1 ## run tests in order defined within this file
random.seed(3141)

# from pandas.testing import assert_frame_equal

class TestGlobalFoldNetworkIrreversible(unittest.TestCase):

    maxDiff = None ## allows full output of failed test differences

    def setUp(self):
        warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        reactions = 10
        compounds = 11
        rids = ['R' + str(x) for x in range(reactions)]
        cids = ['C' + str(x) for x in range(compounds)]
        folds = ['F' + str(x) for x in range(reactions)]
        network = {'rn':[],'direction':[],'cid':[],'s':[]}
        fold_rules  = {'rn': [],'rule':[]}
        i = 0
        for r in rids:
            network['rn'].append(r)
            network['direction'].append('forward')
            network['cid'].append(cids[i])    
            network['s'].append(-1)
            
            network['rn'].append(r)
            network['direction'].append('forward')
            network['cid'].append(cids[i+1])    
            network['s'].append(1)
            fold_rules['rn'].append(r)
            fold_rules['rule'].append(folds[i])
            
            i = i +1

        self.network = pd.DataFrame(network)
        self.fold_rules = pd.DataFrame(fold_rules)
        self.fold_rules['fold_sets'] = self.fold_rules.rule.apply(lambda x: set(x.split('_')))
        self.rn2rules = {d["rn"]:{frozenset(d["fold_sets"])} for d in self.fold_rules.to_dict(orient="records")}

        ## Create Metabolism
        self.met = ne.GlobalMetabolicNetwork(metabolism="dev")
        self.met.network = self.network

        # print(self.met.network)

    def test_GlobalFoldNetwork_create_foldrules2rn(self):
        # warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        expected_ruleids = set([
            ("R0", frozenset({'F0'})),
            ("R1", frozenset({'F1'})),
            ("R2", frozenset({'F2'})),
            ("R3", frozenset({'F3'})),
            ("R4", frozenset({'F4'})),
            ("R5", frozenset({'F5'})),
            ("R6", frozenset({'F6'})),
            ("R7", frozenset({'F7'})),
            ("R8", frozenset({'F8'})),
            ("R9", frozenset({'F9'}))
            ])

        self.assertEqual(set(foldrules.ids), expected_ruleids)

    def test_sort_remaining_foldsets_by_size(self):
        # warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set([]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        expected_foldsets_by_size = {
            1: [frozenset({'F0'}),
                frozenset({'F1'}),
                frozenset({'F2'}),
                frozenset({'F3'}),
                frozenset({'F4'}),
                frozenset({'F5'}),
                frozenset({'F6'}),
                frozenset({'F7'}),
                frozenset({'F8'}),
                frozenset({'F9'})]}

        self.assertEqual(expected_foldsets_by_size, fm.sort_remaining_foldsets_by_size(fm.seed.folds))

    def test_sort_foldsets_by_size(self):
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set([]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)

        foldsets = [
            frozenset({'F6','F1'}),
            frozenset({'F0'}),
            frozenset({'F5'}),
            frozenset({'F1', 'F7'}),
            frozenset({'F3'}),
            frozenset({'F7'}),
            frozenset({'F8'}),
            frozenset({'F9'}),
            frozenset({'F4'}),
            frozenset({'F2'})]

        expected_foldsets_by_size = {
            1: [frozenset({'F0'}),
                frozenset({'F2'}),
                frozenset({'F3'}),
                frozenset({'F4'}),
                frozenset({'F5'}),
                frozenset({'F7'}),
                frozenset({'F8'}),
                frozenset({'F9'})],
            2: [frozenset({'F1', 'F6'}),
                frozenset({'F1', 'F7'})]}

        self.assertEqual(expected_foldsets_by_size, fm.sort_foldsets_by_size(foldsets))

    def test_loop_through_remaining_foldsets(self):
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set([]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        current = nf.Params(folds=fm.seed.folds, cpds=fm.seed.cpds, rns=fm.seed.rns, rules=fm.scope.rules.subset_from_folds(fm.seed.folds))
        size2foldsets = fm.sort_remaining_foldsets_by_size(current.folds)
        key_to_maximize = "rules"
        
        max_effects = fm.loop_through_remaining_foldsets(size2foldsets, current, key_to_maximize)

        self.assertTrue(frozenset({"F0"}) in max_effects)
        self.assertEqual(max_effects[frozenset({'F0'})].cpds, {'C1', 'C0'})
        self.assertEqual(max_effects[frozenset({'F0'})].rns, {'R0'})
        self.assertEqual(max_effects[frozenset({'F0'})].folds, set(['F0']))
        self.assertEqual(max_effects[frozenset({'F0'})].rules.ids, {('R0', frozenset({'F0'}))} )

    def test_select_next_foldset(self):
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set([]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        current = nf.Params(folds=fm.seed.folds, cpds=fm.seed.cpds, rns=fm.seed.rns, rules=fm.scope.rules.subset_from_folds(fm.seed.folds))
        size2foldsets = fm.sort_remaining_foldsets_by_size(current.folds)
        algorithm = "look_ahead_rules"

        next_foldset, effects = fm.select_next_foldset(algorithm, size2foldsets, current)
        self.assertEqual(next_foldset, frozenset({'F0'}))

    def test_FoldMetabolism_rule_order_C0_no_indepdendent(self):
        warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set([]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        result = fm.rule_order(algorithm="look_ahead_rules")

        expected_cpds = {'C0': 0,
                        'C1': 2,
                        'C2': 3,
                        'C3': 4,
                        'C4': 5,
                        'C5': 6,
                        'C6': 7,
                        'C7': 8,
                        'C8': 9,
                        'C9': 10,
                        'C10': 11}
        expected_rns = {'R0': 2,
                        'R1': 3,
                        'R2': 4,
                        'R3': 5,
                        'R4': 6,
                        'R5': 7,
                        'R6': 8,
                        'R7': 9,
                        'R8': 10,
                        'R9': 11}
        expected_folds = {'fold_independent': 0,
                        'F0': 2,
                        'F1': 3,
                        'F2': 4,
                        'F3': 5,
                        'F4': 6,
                        'F5': 7,
                        'F6': 8,
                        'F7': 9,
                        'F8': 10,
                        'F9': 11}

        expected_rules = {("R0", frozenset({'F0'})):2,
                        ("R1", frozenset({'F1'})):3,
                        ("R2", frozenset({'F2'})):4,
                        ("R3", frozenset({'F3'})):5,
                        ("R4", frozenset({'F4'})):6,
                        ("R5", frozenset({'F5'})):7,
                        ("R6", frozenset({'F6'})):8,
                        ("R7", frozenset({'F7'})):9,
                        ("R8", frozenset({'F8'})):10,
                        ("R9", frozenset({'F9'})):11}

        self.assertEqual(expected_cpds, result.cpds)
        self.assertEqual(expected_rns, result.rns)
        self.assertEqual(expected_folds, result.folds)
        self.assertEqual(expected_rules, result.rules)
        
    def test_FoldMetabolism_rule_order_C0_independent_R0R1(self):
        random.seed(3141)
        warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set(["R0","R1"]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        result = fm.rule_order(algorithm="look_ahead_rules", debug=False, ordered_outcome=True)

        expected_cpds = {'C0': 0,
                        'C1': 1,
                        'C2': 1,
                        'C3': 4,
                        'C4': 5,
                        'C5': 6,
                        'C6': 7,
                        'C7': 8,
                        'C8': 9,
                        'C9': 10,
                        'C10': 11}
        
        expected_rns = {'R0': 0,
                        'R1': 0,
                        'R2': 4,
                        'R3': 5,
                        'R4': 6,
                        'R5': 7,
                        'R6': 8,
                        'R7': 9,
                        'R8': 10,
                        'R9': 11}

        expected_folds = {'fold_independent': 0, 
                        'F0': 2,  
                        'F1': 3,  
                        'F2': 4,
                        'F3': 5,
                        'F4': 6,
                        'F5': 7,
                        'F6': 8,
                        'F7': 9,
                        'F8': 10,
                        'F9': 11}

        expected_rules = {("R0", frozenset({'F0'})):2,
                        ("R1", frozenset({'F1'})):3,
                        ("R2", frozenset({'F2'})):4,
                        ("R3", frozenset({'F3'})):5,
                        ("R4", frozenset({'F4'})):6,
                        ("R5", frozenset({'F5'})):7,
                        ("R6", frozenset({'F6'})):8,
                        ("R7", frozenset({'F7'})):9,
                        ("R8", frozenset({'F8'})):10,
                        ("R9", frozenset({'F9'})):11}

        # self.assertEqual(expected_cpds, result.cpds)
        self.assertEqual(expected_rns, result.rns)
        self.assertEqual(expected_folds, result.folds)
        self.assertEqual(expected_rules, result.rules)

    def test_write_results(self):

        random.seed(3141)
        warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)
        foldrules = nf.FoldRules.from_rn2rules(self.rn2rules)
        seed = nf.Params(
            rns = set(["R0","R1"]),
            cpds = set(['C0']),
            folds = set([])
        )

        fm = nf.FoldMetabolism(self.met, foldrules, seed)
        result = fm.rule_order(algorithm="look_ahead_rules", write=True, ordered_outcome=True)

        expected_cpds = {'C0': 0,
                        'C1': 1,
                        'C2': 1,
                        'C3': 4,
                        'C4': 5,
                        'C5': 6,
                        'C6': 7,
                        'C7': 8,
                        'C8': 9,
                        'C9': 10,
                        'C10': 11}
        
        expected_rns = {'R0': 0,
                        'R1': 0,
                        'R2': 4,
                        'R3': 5,
                        'R4': 6,
                        'R5': 7,
                        'R6': 8,
                        'R7': 9,
                        'R8': 10,
                        'R9': 11}

        expected_folds = {'fold_independent': 0, 
                        'F0': 2,  
                        'F1': 3,  
                        'F2': 4,
                        'F3': 5,
                        'F4': 6,
                        'F5': 7,
                        'F6': 8,
                        'F7': 9,
                        'F8': 10,
                        'F9': 11}

        expected_rules = {("R0", frozenset({'F0'})):2,
                        ("R1", frozenset({'F1'})):3,
                        ("R2", frozenset({'F2'})):4,
                        ("R3", frozenset({'F3'})):5,
                        ("R4", frozenset({'F4'})):6,
                        ("R5", frozenset({'F5'})):7,
                        ("R6", frozenset({'F6'})):8,
                        ("R7", frozenset({'F7'})):9,
                        ("R8", frozenset({'F8'})):10,
                        ("R9", frozenset({'F9'})):11}

        final_result = pd.read_pickle(result.final_path)
        temp_result = pd.read_pickle(result.temp_path)

        self.assertEqual(final_result.cpds, expected_cpds)
        self.assertEqual(final_result.rns, expected_rns)
        self.assertEqual(final_result.folds, expected_folds)
        self.assertEqual(final_result.rules, expected_rules)

        self.assertEqual(final_result.cpds, temp_result.cpds)
        self.assertEqual(final_result.rns, temp_result.rns)
        self.assertEqual(final_result.folds, temp_result.folds)
        self.assertEqual(final_result.rules, temp_result.rules)

#### Tests to add
# - test ordering of rules when there's an overlap with seed folds

#########################################################################
## Below tests not updated recently

#     def test_FoldMetabolism_rule_order_C0_independent_R3R5(self):
#         fold_independent_rns = set(["R3","R5"])
#         foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldrules)
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
#         foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldrules)
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
#         foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldrules)
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
#         foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldrules)
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
#         foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
#         fm = nf.FoldMetabolism(self.met, foldrules)
#         fm.seed_cpds = set(['C0'])
#         fm.seed_folds = set([])

#     def test_exit_at_correct_iteration_1(self):
#         pass