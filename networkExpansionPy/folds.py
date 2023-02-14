import pandas as pd
from pathlib import PurePath, Path
from copy import copy, deepcopy
import timeit
from pprint import pprint
from collections import Counter
import random
import itertools

asset_path = PurePath(__file__).parent / "assets"

def rule2rn(rn2rule):
    """
    Returns dictionary mapping of rule:rns from a dictionary of rn:rules

    :param rn2rule: dict mapping rns to rules
    :return: dict mapping rules to rns
    """
    rule2rn = dict()
    for rn, rules in rn2rule.items():
        for fs in rules:
            if fs in rule2rn:
                rule2rn[fs].add(rn)
            else:
                rule2rn[fs] = set([rn])
    return rule2rn

def subset_rule2rn_from_folds(folds, rule2rn):
    """
    Returns a dictionary of rules:rns enabled by folds
    
    :param folds: collection of folds used to subset the rules dict
    :param rules: dict of rule:rns mappings to subset from
    :return: dictionary of rules:rns enabled by folds
    """
    return {k:v for k,v in rule2rn.items() if k <= set(folds)}

def update_iteration_dict(iteration_dict, current, iteration):
    """ 
    Adds the iteration at which new compounds, reactions, and folds appear
        throughout the fold/rule expansion.

    :param iteration_dict: existing dictionary containing keys of "cpds", "rns", and "folds"
                           with values of dictionaries mapping cids, rids, or fids to iteration
                           e.g. {"cpds":{"C00343":0,...},...}
    :param current: a dictionary keyed with "cpds", "rns", and "folds", whose values are sets of 
                    ids that the expansion has reached
    :param iteration: integer of current expansion iteration
    :return: updated iteration_dict
    """
    for dtype, ids in current.items():
        for i in ids:
            if i not in iteration_dict[dtype]:
                iteration_dict[dtype][i] = iteration
    return iteration_dict

########################################################################################################################
########################################################################################################################
# class FoldRules:
#     """
#     Stores how fold rules and reactions map to each other, and which reactions are independent of folds.

#     Invariant to seed compounds, seed folds, etc. This is the universe of possible fold rules. 
#     """

#     def __init__(self, rn2rules, fold_independent_rns):
#         """
#         Construct fold rules from an `rn2rules` dictionary, and a set of `fold_independent_rns`.

#         :param rn2rules: dict of rn:{rules}, where rules are frozensets of folds
#         :param fold_independent_rns: a set of reactions which can occur without folds
#         """

#         self.rn2rules = rn2rules ## all rns to fold rules
#         self.rule2rns = rule2rn(rn2rules) ## all fold rules to rns
#         self.rns = set(rn2rules.keys()) ## reactions in fold network only
#         self.folds = set([i for fs in self.rule2rns.keys() for i in fs]) ## all folds
#         self.fold_independent_rns = fold_independent_rns

class Rule:

    id_iter = itertools.count()

    def __init__(self,rn,foldset):
        self.id = next(self.id_iter)
        self.rn = rn
        self.foldset = foldset

    def __repr__(self):
        return "id:\t\t{0} \nrn:\t\t{1} \nfoldset:\t{2}".format(self.id, self.rn, self.foldset)

class FoldRules:

    def __init__(self,rules:list):
        self.rules = rules
        self.rn2rule = self.rn2rule()
        self.fs2rule = self.fs2rule()
        self.folds = self.folds()

    def __repr__(self):
        return "[\n"+",\n".join([str(i) for i in self.rules])+"]"

    @classmethod
    def from_rn2rules(cls, rn2rules):
        rules = list()
        for rn, fs_set in rn2rules.items():
            for fs in fs_set:
                rules.append(Rule(rn, fs))

        return cls(rules)

    def rns(self):
        return set([r.rn for r in self.rules])

    def folds(self):
        return set([f for r in self.rules for f in r.foldset])

    def foldsets(self):
        return set()

    def subset_from_rns(self, rns):
        # return Rules([r for r in self.rules if r.rn in rns])
        return Rules([self.rn2rule[r] for r in rns])

    def subset_from_folds(self, folds):
        return Rules([r for r in self.rules if r.foldset<=folds])
        # return Rules([self.fs2rule[r] for r in rns])

    def rn2rule(self):
        return {r.rn:r for r in self.rules}

    def fs2rule(self):
        return {r.foldset:r for r in self.rules}

class FoldMetabolism:
    ## Rename to something like FoldScope?
    """
    A class to do fold expansion from a metabolism, foldrules, and seeds.
    """

    def __init__(self, metabolism, foldrules, fold_independent_rns, preexpansion=False):
        """
        Calculates expansion scope after seed compounds are defined. 

        :param metabolism: a GlobalMetabolicNetwork object that defines the compound/reaction rules of network expansion
        :param foldrules: a FoldRules object that defines the fold rules of fold expansion
        :kwarg preexpansion: NOT YET IMPLEMENTED -- but the idea is to allow n steps of regular network expansion before beginning fold expansion
        """
        
        self._m = metabolism ## GlobalMetabolicNetwork object
        self._f = foldrules # FoldRules object

        self.fold_independent_rns = None
        self.seed_folds = None
        self._seed_cpds = None

        self.scope_cpds = None
        self.scope_rns = None
        self.scope_folds = None
        self.scope_rules2rn = None
        self.scope_rn2rules = None

    ## Disallow changing metabolism or foldrules after initialization b/c no setter
    @property 
    def metabolism(self):
        return self._m
    
    @property 
    def m(self):
        return self._m
    
    @property 
    def foldrules(self):
        return self._f
    
    @property 
    def f(self):
        return self._f

    ## Changing the seed_cpds after initializations recalcuates scope properties
    @property
    def seed_cpds(self):
        return self._seed_cpds

    @seed_cpds.setter
    def seed_cpds(self, seed_cpds):
        """
        Calculates properties associated with seed compounds if setting to a new set of seeds
        """
        if (self._seed_cpds == None) or (self._seed_cpds != seed_cpds):
            self._seed_cpds = seed_cpds 
            self.scope_cpds, self.scope_rns = self.calculate_scope(self._seed_cpds)
            self.scope_rules = self.f.rules.subset_from_rns(self.scope_rns)
            # self.scope_rn2rules = {k:v for k,v in self.f.rn2rules.items() if k in self.scope_rns}
            # self.scope_rules2rn = rule2rn(self.scope_rn2rules)
            self.scope_folds = set([i for fs in self.scope_rules2rn.keys() for i in fs])

        else:
            pass 
                
    def calculate_scope(self, seed_cpds):
        """
        Calculate the scope of the seed_cpds with all reactions enabled by the global fold network (including fold_independent reactions)

        :param seed_cpds: collection of compound ids
        :return: set of compounds and set of reactions in the scope
        """
        rn_tup_set = set(self.m.rxns2tuple((self.f.rns | self.f.fold_independent_rns)))
        scope_cpds, scope_rns = self.m.expand(seed_cpds, reaction_mask=rn_tup_set)
        return set(scope_cpds), set([i[0] for i in scope_rns])

    def fold_expand(self, rule2rns, cpds):
        """
        Returns a set of compounds and set of reactions enabled by expansion 
        using `cpds` and the reactions from `rule2rns` and fold_independent 
        reactions.

        :param cpds: collection of compound ids
        :return: set of compounds and set of reactions from the expansion
        """
        fold_rns = set([i for s in rule2rns.values() for i in s])
        rn_tup_set = set(self.m.rxns2tuple(fold_rns | self.f.fold_independent_rns))
        cx,rx = self.m.expand(cpds, reaction_mask=rn_tup_set)
        return set(cx), set([i[0] for i in rx])

    def effect_per_rule_or_fold(self, rule, current_folds, current_cpds):
        """
        Returns the outcomes of an expansion using current folds plus
            folds from an additional `rule`.

        :rule: collection of folds representing the rule
        :current_folds: collection of current folds
        :current_cpds: collection of current compounds
        :return: dict of rule2rns, a set of compounds, and a set of reactions from 
                 running an expansion using the current folds and the `rule` folds. 
        """

        potential_fold_set = (current_folds | set(rule))
        potential_rule2rns = subset_rule2rn_from_folds(potential_fold_set, self.scope_rules2rn)
        cx,rx = self.fold_expand(potential_rule2rns, current_cpds)
        return potential_rule2rns, cx, rx

    ######

    def fold_expand2(self, fold_rns, current_cpds):

        rn_tup_set = set(self.m.rxns2tuple(fold_rns | self.fold_independent_rns))
        cx,rx = self.m.expand(current_cpds, reaction_mask=rn_tup_set)
        return set(cx), set([i[0] for i in rx])


    def effect_per_rule_or_fold2(self, rule, current_folds, current_cpds):

        potential_fold_set = (current_folds | set(rule))
        potential_rules = self.f.subset_from_folds(potential_fold_set)
        cx,rx = self.fold_expand2(potential_rules.rns, current_cpds)
        return cx, rx


    def loop_through_rules2(self):

        potential_ruleset = self.effect_per_rule_or_fold

    ######

    def loop_through_rules(self, current_cpds, current_rns, current_folds, remaining_rules):
        ## future_rule2rns excludes rules which don't lead to new reactions
        future_rule2rns = {k:(v | current_rns) for k,v in remaining_rules.items() if len(v-current_rns)>0}
        rule_sizes = sorted(set([len(i) for i in future_rule2rns]))

        ## Organizes rules by size
        future_rule2rns_by_size = {size:dict() for size in rule_sizes}
        for rule, rns in future_rule2rns.items():
            future_rule2rns_by_size[len(rule)][rule]=rns

        n_rules_checked = 0
        n_rules_skipped = 0
        max_r_effects = dict()
        max_v = 0
        rn_sets_enabling_less_than_max = set()
        for size in rule_sizes:
            ## go through each rule, sorted by maximum rule size
            ## this way i can skip rules if the bigger rule produces no rew reactions, or less than max reactions
            future_rule2rns_by_size_sorted = sorted(future_rule2rns_by_size[size], key=lambda k: len(future_rule2rns_by_size[size][k]), reverse=True)
            for rule in future_rule2rns_by_size_sorted:
                rns = future_rule2rns[rule]
                
                ## Shortcut expansion if the rule only maps to rns which are known not to expand to max
                skip_expansion=False
                for bum_rset in rn_sets_enabling_less_than_max:
                    if rns <= bum_rset:
                        skip_expansion=True
                        print("skipping expansion for ", rule)
                        n_rules_skipped+=1
                        break
                
                ## Expansion
                if not skip_expansion:
                    
                    _fdict = dict()
                    _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(rule, current_folds, current_cpds)

                    n_rules_checked+=1
                    n_new_rns = len(_fdict["rns"] - current_rns)
                    print("n_rules_checked: ", n_rules_checked)
                    if n_new_rns > 0:
                        print("rule enabling new reactions found! ", rule)                    
                    
                    if n_new_rns == max_v:
                        max_r_effects[rule] = _fdict
                    elif n_new_rns > max_v:
                        max_v = n_new_rns
                        for rule in max_r_effects:
                            rn_sets_enabling_less_than_max.add(frozenset(max_r_effects[rule]["rns"]))
                        max_r_effects = dict()
                        max_r_effects[rule] = _fdict
                    else: # n_new_rns < max_v
                        rn_sets_enabling_less_than_max.add(frozenset(_fdict["rns"]))

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(max_r_effects)>0:
                break

        return max_r_effects, n_rules_checked, n_rules_skipped

    def select_next_rule_or_fold(self, current_folds, current_cpds, current_rns, remaining_folds, remaining_rules, algorithm="maxreactionsupersets"):#rselect_func=maxreactions):
        """
        Determines the next rule and its effects on cpds, rns, and rules, based on the selection function.

        :param current_folds: collection of current folds
        :param current_cpds: collection of current compounds
        :param current_rns: collection of current reactions
        :param rselect_func: function to use for identifying the next rule
        :return:- next_rule and ...
                - ... dictionary of its effects {rule2rn:... , cpds:... , rns:...} 
                - and for metadata purposes, also return the number of rules checked
                - and the number of equal rule groups this iteration
        """
        if algorithm == "randomfold":
            next_fold = random.choice(remaining_folds)
            _fdict = dict()
            _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(next_fold, current_folds, current_cpds)
            return next_fold, _fdict, 1, 0

        elif algorithm == "randomrule":
            rule_size_dict = rule_sizes(remaining_rules)
            next_rule = random.choice(rule_size_dict[min(rule_size_dict)])
            _fdict = dict()
            _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(next_fold, current_folds, current_cpds)
            return next_rule, _fdict, 1, 0
            
        elif algorithm == "maxreactionsupersets":
            r_effects, n_rules_checked, n_rules_skipped = self.loop_through_rules(current_cpds, current_rns, current_folds, remaining_rules)
            if len(r_effects) == 0:
                next_rule = frozenset()
                r_effects[next_rule] = {"cpds":deepcopy(current_cpds), "rns":deepcopy(current_rns)}
                print("NO R_EFFECTS REMAINING")
            else:
                next_rule = random.choice(sorted(r_effects.keys()))
            return next_rule, r_effects[next_rule], n_rules_checked, n_rules_skipped #, er_effects

        ## Assuming you start with all scope compounds in the seed set, pick fold order based on fold which enables the most reactions
        elif algorithm == "allcompoundsseedmaxreactions":
            # for r in remaining_rules:
            #     potential_fold_set = (current_folds | set(rule))
            #     potential_rule2rns = subset_rule2rn_from_folds(potential_fold_set, self.scope_rules2rn)

            future_rule2rns = {k:(v | current_rns) for k,v in remaining_rules.items() if len(v-current_rns)>0}
            rule_sizes = sorted(set([len(i) for i in future_rule2rns]))

            ## Organizes rules by size
            future_rule2rns_by_size = {size:dict() for size in rule_sizes}
            for rule, rns in future_rule2rns.items():
                future_rule2rns_by_size[len(rule)][rule]=rns

            for size in rule_sizes:
                future_rule2rns_by_size_sorted = sorted(future_rule2rns_by_size[size], key=lambda k: len(future_rule2rns_by_size[size][k]), reverse=True)
                # future_rule2rns[rule] for rule in future_rule2rns_by_size_sorted:
                    # rns = future_rule2rns[rule]


            ## Don't look for longer rules if shorter rules enable new reactions
            # if len(max_r_effects)>0:
            #     break

        ## picks folds maximizing the number of rules
        elif algorithm == "allcompoundsseedmaxrules":

            rules_minus_current_folds = {k-current_folds:v for k,v in remaining_rules.items()}
            rule_sizes = sorted(set([len(i) for i in rules_minus_current_folds]))
            # print(rule_sizes)
            ## Do we still want rule to be added even if it only maps to reactions which are already in the network?
            ## case 1: want to find rule regardless of it mapping to reactions already in the network

            ## case 2: want to find rule only if it maps to reactions not in the network
            # pprint(rules_minus_current_folds)

            ## pick fold that immeasiately enables the most rules 
            fold2activatedRuleCounter = Counter([f for k in rules_minus_current_folds for f in k if len(k)==min(rule_sizes)])
            fold2activatedRuleMaxValue = max(fold2activatedRuleCounter.values())
            fold2activatedRuleMaxFold = [k for k,v in fold2activatedRuleCounter.items() if v==fold2activatedRuleMaxValue]
            print(fold2activatedRuleMaxFold)
            next_rule = random.choice(sorted(r_effects.keys()))

            ## need to make a conversion between fold2rule and rule2fold

        ## Need to choose an option that picks folds which maximize the number of rules
        # elif algorithm == "maxrules":
        #     blah blah blah

        # elif algorithm == "maxreactions":
        #     r_effects, n_rules_checked, n_equal_rule_groups = self.loop_through_rules(current_folds, current_cpds, current_rns)

    class Metadata:

        def __init__(self):
            self.runtime = dict()
            self.n
    

    
    def rule_order(self, algorithm="maxreactionsupersets"):
        """
        Determine the ordering of all rules/folds.

        :kwarg track_free_rules: if True, add rules/folds to the iteration_dict that
                           weren't selected, but whose reactions are all
                           already enabled.

                           (Probably I want to just remove this kwarg and
                            always have this enabled...)
        """
        valid_algorithms = ["randomfold", "randomrule", "maxreactions", "maxreactionsupersets", "allcompoundsseedmaxreactions", "allcompoundsseedmaxrules"]
        if algorithm.lower() not in valid_algorithms:
            raise ValueError("algorithm must be one of %s"%valid_algorithms)

        if (self.seed_cpds == None) or (self.seed_folds == None):
            raise ValueError("self.seed_cpds and self.seed_folds must not be None")
        
        ## Data to store meta information about run
        metadict = {
            "runtime": dict(),
            "freefolds": dict(),
            "n_rules_skipped":dict(),
            "n_rules_checked":dict(),
            "max_n_remaining_folds":dict(),
            "max_n_remaining_rules":dict()
        }

        ## Initialize current values
        current = {
            "folds": deepcopy(self.seed_folds),
            "cpds": deepcopy(self.seed_cpds),
            "rns":set([])
            }
        
        ## Initialize iteration data 
        iteration = 0
        iteration_dict = {
            "cpds":dict(), 
            "rns":dict(), 
            "folds":{"fold_independent":0}
            }
        
        ################################################
        ## ITERATION 0 (Avoid updating folds on the 0th iteration since they don't apply until iteration=1)
        iteration_dict = update_iteration_dict(iteration_dict, {k:v for k,v in current.items() if k!="folds"}, iteration)
        iteration+=1
        start = timeit.default_timer()

        ################################################
        ## ITERATION 1 (using only seed folds and fold independent reactions)
        init_rules2rn = subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn)
        current["cpds"], current["rns"] = self.fold_expand(init_rules2rn, current["cpds"])
        remaining_folds = (self.scope_folds - current["folds"]) 
        remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if len(k & remaining_folds)>0}
        iteration_dict = update_iteration_dict(iteration_dict, current, iteration)
        ## Update metadata
        metadict["runtime"][iteration] = timeit.default_timer() - start
        # metadict["freefolds"][iteration] = free_folds
        metadict["n_rules_skipped"][iteration] = 0
        metadict["n_rules_checked"][iteration] = 0
        metadict["max_n_remaining_folds"][iteration] = len(remaining_folds)
        metadict["max_n_remaining_rules"][iteration] = len(self.scope_rules2rn) - len(subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn))

        ## Needed in case expansion not possible at all
        if len(remaining_folds) > 0:
            keepgoing = True
        else:
            keepgoing = False

        ## Could it be that there could still be remaining_folds here, yet when we go to select the next rule, none of the folds
        ##      add new reactions? But if this was the case, shouldn't all of those folds be discovered via the free_rules() func?

        # print(f'{free_folds = }')
        # print(f'{remaining_folds = }')
        # remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if k not in subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn)}
        # print(f'{remaining_rules = }')

        ################################################
        ## ITERATION 2+
        while keepgoing:
            start = timeit.default_timer()
            iteration += 1
            print("rule_order iteration: ", iteration)
            for k,v in metadict.items():
                print(k, v)
            next_rule, fdata, n_rules_checked, n_rules_skipped = self.select_next_rule_or_fold(current["folds"], current["cpds"], current["rns"], remaining_folds, remaining_rules, algorithm)
            remaining_folds = (remaining_folds - set(next_rule))
            remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if len(k & remaining_folds)>0}
            # remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if k not in subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn)}

            ## Stop conditions
            if len(remaining_folds) == 0:
                keepgoing = False
            if len(next_rule) == 0:
                keepgoing == False
            if (fdata["cpds"] == current["cpds"]) and (fdata["rns"] == current["rns"]):
                keepgoing = False    
            else:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                current["folds"] = (current["folds"] | set(next_rule))
                current["cpds"] = fdata["cpds"]
                current["rns"] = fdata["rns"]
                
                ## Store when cpds and rns appear in the expansion
                iteration_dict = update_iteration_dict(iteration_dict, current, iteration)

            ## Update metadata
            metadict["runtime"][iteration] = timeit.default_timer() - start
            # metadict["freefolds"][iteration] = free_folds
            metadict["n_rules_skipped"][iteration] = n_rules_skipped
            metadict["n_rules_checked"][iteration] = n_rules_checked
            metadict["max_n_remaining_folds"][iteration] = len(remaining_folds)
            metadict["max_n_remaining_rules"][iteration] = len(self.scope_rules2rn) - len(subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn))

        return current, iteration_dict, metadict

def example_main():
    ## Load metabolism
    metabolism_path = PurePath(asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_path)
    
    ## Load fold rules
    # rn2rules_db_path = PurePath("data", "rn2fold", "_db.pkl")
    # rn2rules_db = pd.read_pickle(rn2rules_db_path)

    # rn2rules_path = rn2rules_db.iloc[4]["OUTPUT_PATH"]
    # rn2rules = pd.read_pickle(rn2rules_path)
    rn2rules_path = PurePath("data", "rn2fold", "rn2rules_v.pkl")
    rn2rules = pd.read_pickle(rn2rules_path)
    foldrules = nf.FoldRules(rn2rules)

    fold_independent_rns = set(metabolism.network["rn"]) - foldrules.rns()
    ## Inititalize fold metabolism
    fm = nf.FoldMetabolism(metabolism, foldrules, fold_independent_rns)

    seed_cpds_path = PurePath("data", "josh", "seed_set.csv")
    fm.seed_cpds = set((pd.read_csv(seed_cpds_path)["ID"]))
    fm.seed_folds = set(['spontaneous'])

    ## Run fold expansion
    current, iteration_dict, metadict = fm.rule_order()
    