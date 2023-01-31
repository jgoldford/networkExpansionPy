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

    # def test_subset_rule2rn(self):
        
    #     rule2rn = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F10','F11'}): {'R0'},
    #         frozenset({'F0'}): {'R0'},
    #         frozenset({'F1'}): {'R1'},
    #         frozenset({'F2'}): {'R2'},
    #         frozenset({'F3'}): {'R3'},
    #         frozenset({'F4'}): {'R4'},
    #         frozenset({'F5'}): {'R5'},
    #         frozenset({'F6'}): {'R6'},
    #         frozenset({'F7'}): {'R7'},
    #         frozenset({'F8'}): {'R8'},
    #         frozenset({'F9'}): {'R9'}
    #         }


    #     folds = {'F10'}
    #     expected = {}
    #     self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

    #     folds = {'F10','F11'}
    #     expected = {
    #         frozenset({'F11'}): {'R2', 'R0'}, 
    #         frozenset({'F10','F11'}): {'R0'}
    #         }
    #     self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

    #     folds = {'F9','F11'}
    #     expected = {
    #         frozenset({'F11'}): {'R2', 'R0'}, 
    #         frozenset({'F9'}): {'R9'}
    #         }
    #     self.assertEqual(expected, nf.subset_rule2rn(folds, rule2rn))

    # def test_rule2nextrns(self):

    #     rule2rn = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F10','F11'}): {'R0'},
    #         frozenset({'F0'}): {'R0'},
    #         frozenset({'F1'}): {'R1'},
    #         frozenset({'F2'}): {'R2'},
    #         frozenset({'F3'}): {'R3'},
    #         frozenset({'F4'}): {'R4'},
    #         frozenset({'F5'}): {'R5'},
    #         frozenset({'F6'}): {'R6'},
    #         frozenset({'F7'}): {'R7'},
    #         frozenset({'F8'}): {'R8'},
    #         frozenset({'F9'}): {'R9'}
    #         }

    #     folds = {'F0'}
    #     expected = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F1'}): {'R1','R0'},
    #         frozenset({'F2'}): {'R2','R0'},
    #         frozenset({'F3'}): {'R3','R0'},
    #         frozenset({'F4'}): {'R4','R0'},
    #         frozenset({'F5'}): {'R5','R0'},
    #         frozenset({'F6'}): {'R6','R0'},
    #         frozenset({'F7'}): {'R7','R0'},
    #         frozenset({'F8'}): {'R8','R0'},
    #         frozenset({'F9'}): {'R9','R0'}
    #         }

    #     self.assertEqual(expected, nf.rule2nextrns(folds, rule2rn))

    #     folds = {'F0', 'F1','F2','F3','F4','F5'}
    #     expected = {
    #         frozenset({'F6'}): {'R6','R0','R1','R2','R3','R4','R5'},
    #         frozenset({'F7'}): {'R7','R0','R1','R2','R3','R4','R5'},
    #         frozenset({'F8'}): {'R8','R0','R1','R2','R3','R4','R5'},
    #         frozenset({'F9'}): {'R9','R0','R1','R2','R3','R4','R5'}
    #         }

    #     self.assertEqual(expected, nf.rule2nextrns(folds, rule2rn))

    # def test_create_equal_rule_groups(self):
    #     rule2rn = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F10','F11'}): {'R0'},
    #         frozenset({'F0'}): {'R0'},
    #         frozenset({'F1'}): {'R1'},
    #         frozenset({'F2'}): {'R2'},
    #         frozenset({'F3'}): {'R3'},
    #         frozenset({'F4'}): {'R4'},
    #         frozenset({'F5'}): {'R5'},
    #         frozenset({'F6'}): {'R6'},
    #         frozenset({'F7'}): {'R7'},
    #         frozenset({'F8'}): {'R8'},
    #         frozenset({'F9'}): {'R9'}
    #         }

    #     expected = {
    #         frozenset([frozenset({'F11'}),frozenset({'F12'})]),
    #         frozenset([frozenset({'F1'})]),
    #         frozenset([frozenset({'F3'})]),
    #         frozenset([frozenset({'F4'})]),
    #         frozenset([frozenset({'F5'})]),
    #         frozenset([frozenset({'F6'})]),
    #         frozenset([frozenset({'F7'})]),
    #         frozenset([frozenset({'F8'})]),
    #         frozenset([frozenset({'F9'})])
    #     }

    #     self.assertEqual(expected, nf.create_equal_rule_groups(rule2rn))

    # def test_sort_equal_rule_groups(self):

    #     equal_rule_groups = {
    #         frozenset([frozenset({'F11'}),frozenset({'F12'})]),
    #         frozenset([frozenset({'F1'})]),
    #         frozenset([frozenset({'F3'})]),
    #         frozenset([frozenset({'F4'})]),
    #         frozenset([frozenset({'F5'})]),
    #         frozenset([frozenset({'F6'})]),
    #         frozenset([frozenset({'F7'})]),
    #         frozenset([frozenset({'F8'})]),
    #         frozenset([frozenset({'F9'})])
    #     }

    #     expected = [
    #         [['F1']],
    #         [['F11'], ['F12']],
    #         [['F3']],
    #         [['F4']],
    #         [['F5']],
    #         [['F6']],
    #         [['F7']],
    #         [['F8']],
    #         [['F9']]
    #     ]

    #     self.assertEqual(expected, nf.sort_equal_rule_groups(equal_rule_groups))

    # def test_rule_sizes(self):
    #     ## LEFT OFF HERE
    #     self.assertEqual({1:[['F1']]}, nf.rule_sizes([['F1']]))
    #     self.assertEqual({2:[['F1','F2']]}, nf.rule_sizes([['F1','F2']]))
    #     self.assertEqual({1:[['F2']], 2:[['F1','F2']]}, nf.rule_sizes([['F1','F2'],['F2']]))

    # def test_remove_current_folds_from_equal_rule_groups(self):

    #     rule_groups = [
    #         [['F1']],
    #         [['F11'], ['F12']],
    #         [['F2','F3']],
    #         [['F4']],
    #         [['F5']],
    #         [['F6']],
    #         [['F7']],
    #         [['F8']],
    #         [['F9']]
    #     ]

    #     current_folds = {'F1'}
    #     expected = [
    #         [frozenset({'F11'}), frozenset({'F12'})],
    #         [frozenset({'F3', 'F2'})],
    #         [frozenset({'F4'})],
    #         [frozenset({'F5'})],
    #         [frozenset({'F6'})],
    #         [frozenset({'F7'})],
    #         [frozenset({'F8'})],
    #         [frozenset({'F9'})]
    #         ]
    #     self.assertEqual(expected, nf.remove_current_folds_from_equal_rule_groups(current_folds, rule_groups))

    #     current_folds = {'F2','F11'}
    #     expected = [
    #         [frozenset({'F1'})],
    #         [frozenset({'F12'})],
    #         [frozenset({'F3'})],
    #         [frozenset({'F4'})],
    #         [frozenset({'F5'})],
    #         [frozenset({'F6'})],
    #         [frozenset({'F7'})],
    #         [frozenset({'F8'})],
    #         [frozenset({'F9'})]
    #         ]
    #     self.assertEqual(expected, nf.remove_current_folds_from_equal_rule_groups(current_folds, rule_groups))

    # def test_next_iter_possible_rules(self):

    #     current_folds = {'F3','F11'}
    #     rule2rn = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F10','F11'}): {'R0'},
    #         frozenset({'F0'}): {'R0'},
    #         frozenset({'F1'}): {'R1'},
    #         frozenset({'F2'}): {'R2'},
    #         frozenset({'F3'}): {'R3'},
    #         frozenset({'F4'}): {'R4'},
    #         frozenset({'F5'}): {'R5'},
    #         frozenset({'F6'}): {'R6'},
    #         frozenset({'F7'}): {'R7'},
    #         frozenset({'F8'}): {'R8'},
    #         frozenset({'F9'}): {'R9'}
    #         }
    #     expected = [
    #         {1:[frozenset({'F1'})]},
    #         {1:[frozenset({'F4'})]},
    #         {1:[frozenset({'F5'})]},
    #         {1:[frozenset({'F6'})]},
    #         {1:[frozenset({'F7'})]},
    #         {1:[frozenset({'F8'})]},
    #         {1:[frozenset({'F9'})]},
    #     ]
    #     self.assertEqual(expected, nf.next_iter_possible_rules(current_folds, rule2rn))

    #     current_folds = {'F10'}
    #     expected = [
    #         {1:[frozenset({'F1'})]},
    #         {1:[frozenset({'F11'}), frozenset({'F12'})]}, ## F11 and F12 enable equivilent rxn sets
    #         {1:[frozenset({'F3'})]},
    #         {1:[frozenset({'F4'})]},
    #         {1:[frozenset({'F5'})]},
    #         {1:[frozenset({'F6'})]},
    #         {1:[frozenset({'F7'})]},
    #         {1:[frozenset({'F8'})]},
    #         {1:[frozenset({'F9'})]},
    #     ]
    #     ## F0 and F2 enable only subsets of F11 and F12 so are excluded
    #     self.assertEqual(expected, nf.next_iter_possible_rules(current_folds, rule2rn))

    # def test_maxreactions(self):
    #     r_effects = {
    #         frozenset({'F1'}):{"rns":[1,2,3,4,5]},
    #         frozenset({'F1','F3'}):{"rns":[1,2,3,4,5,4,1,2,3,4,]},
    #         frozenset({'F90'}):{"rns":[1,2]}
    #     }

    #     self.assertEqual(frozenset({'F1','F3'}), nf.maxreactions(r_effects))

    # def test_update_iteration_dict(self):
    #     ## skip testing for now
    #     pass

    # def test_free_rules(self):
    #     current_folds = {'F3','F11'}
    #     rule2rn = {
    #         frozenset({'F12'}): {'R2', 'R0'},
    #         frozenset({'F11'}): {'R2', 'R0'},
    #         frozenset({'F10','F11'}): {'R0'},
    #         frozenset({'F0'}): {'R0'},
    #         frozenset({'F1'}): {'R1'},
    #         frozenset({'F2'}): {'R2'},
    #         frozenset({'F3'}): {'R3'},
    #         frozenset({'F4'}): {'R4'},
    #         frozenset({'F5'}): {'R5'},
    #         frozenset({'F6'}): {'R6'},
    #         frozenset({'F7'}): {'R7'},
    #         frozenset({'F8'}): {'R8'},
    #         frozenset({'F9'}): {'R9'}
    #         }

    #     expected = {frozenset({'F10','F11'}), frozenset({'F12'}), frozenset({'F0'}), frozenset({'F2'})}

    #     self.assertEqual(expected, nf.free_rules(current_folds, rule2rn))

    # def test_next_iter_possible_rules(self):
    #     current_folds = {'spontaneous', '246', '210', '2007', '286'}
    #     remaining_rules = {frozenset({'7579'}): {'R09126', 'R09983'},
    #         frozenset({'7561'}): {'R03540'},
    #         frozenset({'149'}): {'R11750'},
    #         frozenset({'101', '812'}): {'R03546', 'R10079'},
    #         frozenset({'315'}): {'R06973'},
    #         frozenset({'4952', '4953', '602'}): {'R00490'},
    #         frozenset({'4021'}): {'R02243'},
    #         frozenset({'4952', '602'}): {'R11749'},
    #         frozenset({'7507', '7546'}): {'R00485'},
    #         frozenset({'7523'}): {'R02244'},
    #         frozenset({'2003', '7574'}): {'R00224', 'R00636'},
    #         frozenset({'7552'}): {'R00321'},
    #         frozenset({'278', '604'}): {'LAO_FM'}}
    #     current_rns = {'R00269',
    #         'R00348',
    #         'R01087',
    #         'R02244',
    #         'R07316',
    #         'R08698',
    #         'R11617',
    #         'R12185'}

    #     ## [(i.equal_supersets,i.equal_subsets) for i in equal_rule_dict] probably need to reformat to sets
    #     #  expected_equal_rule_groups = [([frozenset({'315'})], []),
    #         # ([frozenset({'2003', '7574'})], []),
    #         # ([frozenset({'278', '604'})], []),
    #         # ([frozenset({'101', '812'})], []),
    #         # ([frozenset({'4952', '602'})], []),
    #         # ([frozenset({'210'})], [frozenset({'7507', '7546'})]),
    #         # ([frozenset({'7561'})], []),
    #         # ([frozenset({'7552'})], []),
    #         # ([frozenset({'149'})], []),
    #         # ([frozenset({'4021'})], []),
    #         # ([frozenset({'7579'})], []),
    #         # ([frozenset({'4952', '4953', '602'})], [])]

    #     equal_rule_groups = nf.next_iter_possible_rules(current["folds"], remaining_rules, current["rns"])

    # def test_loop_through_rules(self):

    #     current_folds = {'2007', '246', '7523', 'spontaneous'}
    #     current_cpds = {'C00001',
    #         'C00009',
    #         'C00011',
    #         'C00012',
    #         'C00014',
    #         'C00022',
    #         'C00023',
    #         'C00026',
    #         'C00028',
    #         'C00030',
    #         'C00033',
    #         'C00036',
    #         'C00038',
    #         'C00042',
    #         'C00048',
    #         'C00058',
    #         'C00069',
    #         'C00070',
    #         'C00071',
    #         'C00080',
    #         'C00122',
    #         'C00149',
    #         'C00150',
    #         'C00160',
    #         'C00161',
    #         'C00175',
    #         'C00205',
    #         'C00209',
    #         'C00238',
    #         'C00282',
    #         'C00283',
    #         'C00288',
    #         'C00305',
    #         'C00311',
    #         'C00383',
    #         'C00417',
    #         'C00940',
    #         'C01127',
    #         'C01330',
    #         'C01335',
    #         'C01384',
    #         'C01528',
    #         'C01563',
    #         'C01732',
    #         'C02218',
    #         'C02341',
    #         'C02362',
    #         'C06232',
    #         'C14818',
    #         'C14819',
    #         'C17023',
    #         'C19609',
    #         'C19806',
    #         'C20679',
    #         'C22155',
    #         'Z00001',
    #         'Z00002',
    #         'Z00006',
    #         'Z00015',
    #         'Z00020',
    #         'Z00029',
    #         'Z00030',
    #         'Z00033',
    #         'Z00034',
    #         'Z00053',
    #         'Z00054',
    #         'Z00055',
    #         'Z00060',
    #         'Z00062',
    #         'Z00063',
    #         'Z00064',
    #         'Z00067',
    #         'Z00069',
    #         'Z00070'}
    #     current_rns = {'R00269', 'R00348', 'R01087', 'R02244', 'R07316', 'R08698', 'R12185'}
    #     remaining_rules = {frozenset({'7579'}): {'R09126', 'R09983'},
    #         frozenset({'7561'}): {'R03540'},
    #         frozenset({'149'}): {'R11750'},
    #         frozenset({'101', '812'}): {'R03546', 'R10079'},
    #         frozenset({'315'}): {'R06973'},
    #         frozenset({'210'}): {'R00485', 'R11617'},
    #         frozenset({'4952', '4953', '602'}): {'R00490'},
    #         frozenset({'4021'}): {'R02243'},
    #         frozenset({'4952', '602'}): {'R11749'},
    #         frozenset({'7507', '7546'}): {'R00485'},
    #         frozenset({'286'}): {'R02244'},
    #         frozenset({'2003', '7574'}): {'R00224', 'R00636'},
    #         frozenset({'7552'}): {'R00321'},
    #         frozenset({'278', '604'}): {'LAO_FM'}}

    #     expected_r_effects_keys = dict_keys([frozenset({'210'}), frozenset({'7561'}), frozenset({'7552'}), frozenset({'4021'})])
    #     # expected_n_rules_checked = 
    #     ## [(i.equal_supersets,i.equal_subsets) for i in equal_rule_dict] probably need to reformat to sets
    #     expected_equal_rule_dict = [([frozenset({'315'})], []),
    #         ([frozenset({'2003', '7574'})], []),
    #         ([frozenset({'278', '604'})], []),
    #         ([frozenset({'101', '812'})], []),
    #         ([frozenset({'4952', '602'})], []),
    #         ([frozenset({'210'})], [frozenset({'7507', '7546'})]),
    #         ([frozenset({'7561'})], []),
    #         ([frozenset({'7552'})], []),
    #         ([frozenset({'149'})], []),
    #         ([frozenset({'4021'})], []),
    #         ([frozenset({'7579'})], []),
    #         ([frozenset({'4952', '4953', '602'})], [])]
        
    #     ## [(i.equal_supersets,i.equal_subsets) for i in er_effects] probably need to reformat to sets
    #     expected_er_effects = [([frozenset({'315'})], []),
    #         ([frozenset({'210'})], [frozenset({'7507', '7546'})]),
    #         ([frozenset({'7561'})], []),
    #         ([frozenset({'7552'})], []),
    #         ([frozenset({'149'})], []),
    #         ([frozenset({'4021'})], []),
    #         ([frozenset({'7579'})], [])]

    #     r_effects, n_rules_checked, equal_rule_dict, er_effects = fm.loop_through_rules(current["folds"], current["cpds"], current["rns"], remaining_rules)

class TestGlobalFoldNetworkIrreversible(unittest.TestCase):

    maxDiff = None ## allows full output of failed test differences

    def setUp(self):
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

        print(self.met.network)

    def test_GlobalFoldNetwork_create_foldrules2rn(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        expected_rule2rns = {frozenset({'F0'}): {'R0'},
                        frozenset({'F1'}): {'R1'},
                        frozenset({'F2'}): {'R2'},
                        frozenset({'F3'}): {'R3'},
                        frozenset({'F4'}): {'R4'},
                        frozenset({'F5'}): {'R5'},
                        frozenset({'F6'}): {'R6'},
                        frozenset({'F7'}): {'R7'},
                        frozenset({'F8'}): {'R8'},
                        frozenset({'F9'}): {'R9'}}
        self.assertEqual(foldrules.rule2rns, expected_rule2rns)

    def test_effect_per_rule_or_fold(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])

        expected_potential_rule2rns = {frozenset({'F0'}): {'R0'}}
        expected_cx = {'C0','C1'}
        expected_rx = {'R0'}

        potential_rule2rns, cx, rx = fm.effect_per_rule_or_fold(frozenset({'F0'}), set(), {"C0"})
        self.assertEqual(expected_potential_rule2rns, potential_rule2rns)
        self.assertEqual(expected_cx, cx)
        self.assertEqual(expected_rx, rx)

        expected_potential_rule2rns = {
            frozenset({'F0'}): {'R0'},
            frozenset({'F1'}): {'R1'}
            }
        expected_cx = {'C0','C1','C2'}
        expected_rx = {'R0','R1'}

        potential_rule2rns, cx, rx = fm.effect_per_rule_or_fold(frozenset({'F1'}), {'F0'}, {"C0","C1"})
        self.assertEqual(expected_potential_rule2rns, potential_rule2rns)
        self.assertEqual(expected_cx, cx)
        self.assertEqual(expected_rx, rx)

    def test_scope_rn2rules(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
                                frozenset({'F1'}): {'R1'},
                                frozenset({'F2'}): {'R2'},
                                frozenset({'F3'}): {'R3'},
                                frozenset({'F4'}): {'R4'},
                                frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
                                'R1': {frozenset({'F1'})},
                                'R2': {frozenset({'F2'})},
                                'R3': {frozenset({'F3'})},
                                'R4': {frozenset({'F4'})},
                                'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)


    def test_FoldMetabolism_rule_order_C0_no_indepdendent(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C0': 0,
                                'C1': 2,
                                'C2': 3,
                                'C3': 4,
                                'C4': 5,
                                'C5': 6,
                                'C6': 7,
                                'C7': 8,
                                'C8': 9,
                                'C9': 10,
                                'C10': 11},
                                'rns': {'R0': 2,
                                'R1': 3,
                                'R2': 4,
                                'R3': 5,
                                'R4': 6,
                                'R5': 7,
                                'R6': 8,
                                'R7': 9,
                                'R8': 10,
                                'R9': 11},
                                'folds': {'fold_independent': 0,
                                'F0': 2,
                                'F1': 3,
                                'F2': 4,
                                'F3': 5,
                                'F4': 6,
                                'F5': 7,
                                'F6': 8,
                                'F7': 9,
                                'F8': 10,
                                'F9': 11}}

        self.assertEqual(iteration_dict, expected_iteration_dict)
        self.assertEqual(current, expected_current)
        
    def test_FoldMetabolism_rule_order_C0_independent_R0R1(self):
        fold_independent_rns = set(["R0","R1"])
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
                                frozenset({'F1'}): {'R1'},
                                frozenset({'F2'}): {'R2'},
                                frozenset({'F3'}): {'R3'},
                                frozenset({'F4'}): {'R4'},
                                frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
                                'R1': {frozenset({'F1'})},
                                'R2': {frozenset({'F2'})},
                                'R3': {frozenset({'F3'})},
                                'R4': {frozenset({'F4'})},
                                'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C0': 0,
                                'C1': 1,
                                'C2': 1,
                                'C3': 2,
                                'C4': 3,
                                'C5': 4,
                                'C6': 5,
                                'C7': 6,
                                'C8': 7,
                                'C9': 8,
                                'C10': 9},
                                'rns': {'R0': 1,
                                'R1': 1,
                                'R2': 2,
                                'R3': 3,
                                'R4': 4,
                                'R5': 5,
                                'R6': 6,
                                'R7': 7,
                                'R8': 8,
                                'R9': 9},
                                'folds': {'fold_independent': 0,
                                'F0': 2,
                                'F1': 3,
                                'F2': 4,
                                'F3': 5,
                                'F4': 6,
                                'F5': 7,
                                'F6': 8,
                                'F7': 9,
                                'F8': 10,
                                'F9': 11}}

    def test_FoldMetabolism_rule_order_C0_independent_R3R5(self):
        fold_independent_rns = set(["R3","R5"])
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
                                frozenset({'F1'}): {'R1'},
                                frozenset({'F2'}): {'R2'},
                                frozenset({'F3'}): {'R3'},
                                frozenset({'F4'}): {'R4'},
                                frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
                                'R1': {frozenset({'F1'})},
                                'R2': {frozenset({'F2'})},
                                'R3': {frozenset({'F3'})},
                                'R4': {frozenset({'F4'})},
                                'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C0': 0,
                                'C1': 2,
                                'C2': 3,
                                'C3': 4,
                                'C4': 4,
                                'C5': 5,
                                'C6': 5,
                                'C7': 6,
                                'C8': 7,
                                'C9': 8,
                                'C10': 9},
                                'rns': {'R0': 2,
                                'R1': 3,
                                'R2': 4,
                                'R3': 4,
                                'R4': 5,
                                'R5': 5,
                                'R6': 6,
                                'R7': 7,
                                'R8': 8,
                                'R9': 9},
                                'folds': {'fold_independent': 0,
                                'F0': 2,
                                'F1': 3,
                                'F2': 4,
                                'F3': 5,
                                'F4': 6,
                                'F5': 7,
                                'F6': 8,
                                'F7': 9,
                                'F8': 10,
                                'F9': 11}}

    def test_FoldMetabolism_rule_order_C5_no_independent(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C5'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C10', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C5': 0, 'C6': 2, 'C7': 3, 'C8': 4, 'C9': 5, 'C10': 6},
                            'rns': {'R5': 2, 'R6': 3, 'R7': 4, 'R8': 5, 'R9': 6},
                            'folds': {'fold_independent': 0, 'F5': 2, 'F6': 3, 'F7': 4, 'F8': 5, 'F9': 6}}
        self.assertEqual(iteration_dict, expected_iteration_dict)
        self.assertEqual(current, expected_current)

class TestGlobalFoldNetworkReversible(unittest.TestCase):

    def setUp(self):
        reactions = 10
        compounds = 11
        rids = ['R' + str(x) for x in range(reactions)]
        cids = ['C' + str(x) for x in range(compounds)]
        folds = ['F' + str(x) for x in range(reactions)]
        network = {'rn':[],'direction':[],'cid':[],'s':[]}
        fold_rules  = {'rn': [],'rule':[]}
        i = 0
        for r in rids:
            ## Forward i
            network['rn'].append(r)
            network['direction'].append('forward')
            network['cid'].append(cids[i])    
            network['s'].append(-1)
            ## Backwards i
            network['rn'].append(r)
            network['direction'].append('backward')
            network['cid'].append(cids[i])    
            network['s'].append(1)
            
            ## Forward i+1
            network['rn'].append(r)
            network['direction'].append('forward')
            network['cid'].append(cids[i+1])    
            network['s'].append(1)
            fold_rules['rn'].append(r)
            fold_rules['rule'].append(folds[i])
            ## Backwards i+1
            network['rn'].append(r)
            network['direction'].append('backward')
            network['cid'].append(cids[i+1])    
            network['s'].append(-1)
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

    def test_FoldMetabolism_rule_order_C0_reversible(self):
        fold_independent_rns = set()
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
                                frozenset({'F1'}): {'R1'},
                                frozenset({'F2'}): {'R2'},
                                frozenset({'F3'}): {'R3'},
                                frozenset({'F4'}): {'R4'},
                                frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
                                'R1': {frozenset({'F1'})},
                                'R2': {frozenset({'F2'})},
                                'R3': {frozenset({'F3'})},
                                'R4': {frozenset({'F4'})},
                                'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C0': 0,
                                'C1': 2,
                                'C2': 3,
                                'C3': 4,
                                'C4': 5,
                                'C5': 6,
                                'C6': 7,
                                'C7': 8,
                                'C8': 9,
                                'C9': 10,
                                'C10': 11},
                                'rns': {'R0': 2,
                                'R1': 3,
                                'R2': 4,
                                'R3': 5,
                                'R4': 6,
                                'R5': 7,
                                'R6': 8,
                                'R7': 9,
                                'R8': 10,
                                'R9': 11},
                                'folds': {'fold_independent': 0,
                                'F0': 2,
                                'F1': 3,
                                'F2': 4,
                                'F3': 5,
                                'F4': 6,
                                'F5': 7,
                                'F6': 8,
                                'F7': 9,
                                'F8': 10,
                                'F9': 11}}

        self.assertEqual(iteration_dict, expected_iteration_dict)
        self.assertEqual(current, expected_current)

    def test_FoldMetabolism_rule_order_C5_reversible(self):
        fold_independent_rns = set([])
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C5'])
        fm.seed_folds = set([])

        expected_scope_rules2rn = {frozenset({'F0'}): {'R0'},
                                frozenset({'F1'}): {'R1'},
                                frozenset({'F2'}): {'R2'},
                                frozenset({'F3'}): {'R3'},
                                frozenset({'F4'}): {'R4'},
                                frozenset({'F5'}): {'R5'},
                                frozenset({'F6'}): {'R6'},
                                frozenset({'F7'}): {'R7'},
                                frozenset({'F8'}): {'R8'},
                                frozenset({'F9'}): {'R9'}}
        self.assertEqual(fm.scope_rules2rn, expected_scope_rules2rn)

        expected_scope_rn2rules = {'R0': {frozenset({'F0'})},
                                'R1': {frozenset({'F1'})},
                                'R2': {frozenset({'F2'})},
                                'R3': {frozenset({'F3'})},
                                'R4': {frozenset({'F4'})},
                                'R5': {frozenset({'F5'})},
                                'R6': {frozenset({'F6'})},
                                'R7': {frozenset({'F7'})},
                                'R8': {frozenset({'F8'})},
                                'R9': {frozenset({'F9'})}}
        self.assertEqual(fm.scope_rn2rules, expected_scope_rn2rules)

        current, iteration_dict, metadict = fm.rule_order()
        expected_current = {'folds': {'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9'},
                            'cpds': {'C0', 'C1', 'C10', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'},
                            'rns': {'R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9'}}

        expected_iteration_dict = {'cpds': {'C0': 6,
                                'C1': 5,
                                'C2': 4,
                                'C3': 3,
                                'C4': 2,
                                'C5': 0,
                                'C6': 7,
                                'C7': 8,
                                'C8': 9,
                                'C9': 10,
                                'C10': 11},
                                'rns': {'R0': 6,
                                'R1': 5,
                                'R2': 4,
                                'R3': 3,
                                'R4': 2,
                                'R5': 7,
                                'R6': 8,
                                'R7': 9,
                                'R8': 10,
                                'R9': 11},
                                'folds': {'fold_independent': 0,
                                'F0': 2,
                                'F1': 3,
                                'F2': 4,
                                'F3': 5,
                                'F4': 6,
                                'F5': 7,
                                'F6': 8,
                                'F7': 9,
                                'F8': 10,
                                'F9': 11}}

class TestGlobalFoldNetworkTwoFoldsSimultaneouslyNeeded(unittest.TestCase):

    def setUp(self):
        reactions = 10
        compounds = 11
        rids = ['R' + str(x) for x in range(reactions)]
        cids = ['C' + str(x) for x in range(compounds)]
        folds = ['F' + str(x) for x in range(reactions)]
        network = {'rn':[],'direction':[],'cid':[],'s':[]}
        fold_rules  = {'rn': [],'rule':[]}
        i = 0
        for r in rids:
            ## Forward i
            network['rn'].append(r)
            network['direction'].append('forward')
            network['cid'].append(cids[i])    
            network['s'].append(-1)

            ## Forward i+1
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
        self.rn2rules["R0"] = {frozenset({'F0','F10'})}

        ## Create Metabolism
        self.met = ne.GlobalMetabolicNetwork(metabolism="dev")
        self.met.network = self.network

    def test_FoldMetabolism_rule_order_R0_needs_2_folds(self):
        fold_independent_rns = set([])
        foldrules = nf.FoldRules(self.rn2rules, fold_independent_rns)
        fm = nf.FoldMetabolism(self.met, foldrules)
        fm.seed_cpds = set(['C0'])
        fm.seed_folds = set([])

    def test_exit_at_correct_iteration_1(self):
        pass