import pandas as pd
from pathlib import PurePath, Path
from copy import deepcopy
import timeit
from pprint import pprint
from collections import Counter
import random

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

# def rule2nextrns(current_folds, rule2rn, restrict_to_new = True):

#     """
#     Returns a dictionary of rules:rns only for rules enabling reactions that 
#         are undiscovered by current_folds. Each rule will map to all reactions
#         enabled after that rule is added.

#     :param current_folds: collection of folds to compare to
#     :param scope_rule2rn: dict of rule:rns mappings to subset from, and then compare to
#     :kwarg restrict_to_new: if true, only return rules which enable undiscoverd reactions
#                             if false, include all rules, even those which don't add any new reactions and which have been discovered
#     """
#     current_rule2rn = subset_rule2rn_from_folds(current_folds, scope_rule2rn)
#     return {k:v for k,v in scope_rule2rn.items() if not k in current_rule2rn} 
#     ## This would only return new rules, regardless of if they contribute new reactions. seeems like the best way to do this.
#     ##      but right now this only maps rules:that_rules_rns, not to all reactions which would be possible in total if that rule was added


#     permitted_rns = set([rn for v in current_rule2rn.values() for rn in v])
#     # print("interest- ", rule2rn[frozenset({'304', '222', '7581', '3321', '2002', '3323'})])
#     # print(frozenset({'304', '222', '7581', '3321', '2002', '3323'}) <= permitted_rns)
#     # for k,v in rule2rn.items():
#     #     if k == frozenset({'304', '222', '7581', '3321', '2002', '3323'}):
#     #         print("found")
#     #         print(k, v)
#     #         print(v <= permitted_rns)
#     #     if not (v <= permitted_rns):
#     #         print(k, v)
#     if restrict_to_new:
#         return {k:(v | permitted_rns) for k,v in scope_rule2rn.items() if not v <= permitted_rns} ## will exclude already discovered rules
#     else:
#         return {k:(v | permitted_rns) for k,v in scope_rule2rn.items()} ## won't exclude already discovered rules

# def create_equal_rule_groups(rule2rn):
#     """
#     Returns a set of equivilent rule collections.

#     Collections and the rules they contain are each frozensets.

#     :param rule2rn: dict of rule:rns to create equal rule groups from
#     :return: set of frozensets (i.e. groups of equivilent rules) of frozensets (i.e. rules) 
#     """

#     strictsubset_ks = set()
#     equal_groups = set()
#     for k, v in rule2rn.items():
#         equal_ks = {k}
#         for k2, v2 in rule2rn.items():
#             if k != k2:
#                 if v2 <= v:
#                     if v2!=v:
#                         strictsubset_ks.add(k2)
#                     else: # v2==v
#                         equal_ks.add(k2)
        
#         equal_groups.add(frozenset(equal_ks))

#     ## exclude groups which include any strict subset rules
#     return {i for i in equal_groups if not i & strictsubset_ks}

def create_equal_rule_groups(rule2rn):
    """
    Returns a set of equivilent rule collections.

    Collections and the rules they contain are each frozensets.

    :param rule2rn: dict of rule:rns to create equal rule groups from
    :return: set of frozensets (i.e. groups of equivilent rules) of frozensets (i.e. rules) 
    """

    strictsubsets = set()
    # equal_groups = set()
    equal_groups = []
    for k, v in rule2rn.items():
        equal_ks = {k}
        strictsubset_ks = set()
        for k2, v2 in rule2rn.items():
            if k != k2:
                if v2 <= v:
                    if v2!=v:
                        # maps to the a strict subset of reactions
                        strictsubsets.add(k2)
                        strictsubset_ks.add(k2)
                    else: # v2==v
                        # maps to the exact same reactions
                        equal_ks.add(k2) 
        
        er = EquivilentRule(equal_ks, strictsubset_ks)
        if er not in equal_groups:
            equal_groups.append(er)
        # equal_groups.append((frozenset(equal_ks),frozenset(strictsubset_ks)))
        # equal_groups.append({frozenset(equal_ks):frozenset(strictsubset_ks)})
    
    # equal_groups = set(equal_groups)
    ## exclude groups which include any strict subset rules
    # return {d for d in equal_groups if not d[0] & strictsubsets}
    # return [{d[0]:d[1]} for d in equal_groups if not d[0] & strictsubsets]
    # return [{"superset":d[0], "subset":d[1]} for d in equal_groups if not d[0] & strictsubsets]
    return set([i for i in equal_groups if not i.equal_supersets & strictsubsets]) ## ignore equivilent rules that have supersets

def sort_equal_rule_groups(equal_rule_groups):
    """
    Returns a sorted list of equivilent rule collections.

    Collections and list are each sorted.

    :param equal_rule_groups: set of equivilent rule collections
    :return: a sorted list of equivilent rule collections (each sorted)
    """
    equal_rule_groups_sorted = []
    for ofs in equal_rule_groups:
        equal_rule_groups_sorted.append(sorted([sorted(ifs) for ifs in ofs]))

    return sorted(equal_rule_groups_sorted)

def rule_sizes(rule_group):
    """
    Returns a dictionary of len:collections for collections in rule_group

    :param rule_group: a set of collections
    :return: dictionary keyed by lengths of each collection in value
    """
    element_lengths = dict()
    for i in rule_group:
        l = len(i)
        if not l in element_lengths:
            element_lengths[l] = [i]
        else:
            element_lengths[l].append(i)
    return element_lengths

# def remove_current_folds_from_equal_rule_groups(current_folds, rule_groups):
#     """ 
#     Returns a list of rule groups with current folds removed.

#     Preserves ordering from `rule_groups`.

#     :param current_folds: collection of current folds
#     :param rule_groups: a collection of equivilent rule groups
#     :return: a list of rule groups with current folds removed
#     """

#     new_rule_groups = []
#     for rg in rule_groups:
#         new_group = []
#         # this construction allows us to preserve the ordering from create_equal_rule_groups
#         for i in rg:
#             new_rule = frozenset(i) - current_folds
#             if len(new_rule) > 0:
#                 new_group.append(new_rule)
#         if len(new_group) > 0:
#             new_rule_groups.append(new_group)

#     return new_rule_groups

# def next_iter_possible_rules(current_folds2, scope_rule2rn, remaining_rules, current_rns):
#     """
#     Returns a list of equal rule group dictionaries, keyed by rule size.

#     :param current_folds: collection of current folds
#     :param scope_rule2rn: dict of rule:rns (should be for scope)
#     :return: a list of equal rule group dictionaries, keyed by rule size
#             [
#                 {1:[rule1,rule2,rule2],3:[rule4,rule5]}, #dict of equal rules, 3 with len=1, 2 with len=3
#                 {1:[rule7,rule8],3:[rule9]}              
#             ]
#     """

#     ## Need to run these two calls every iteration of the fold expansion
#     # future_rule2rns = rule2nextrns(current_folds, scope_rule2rn) ## This is just finnding rules from reamaining folds with current reactions appended to all values
#     future_rule2rns = {k:(v | current_rns) for k,v in remaining_rules.items()}
#     equal_rule_groups = create_equal_rule_groups(future_rule2rns) ## discards any rules which are strict subsets of others...
#     equal_rule_groups = sort_equal_rule_groups(equal_rule_groups)
#     equal_rule_groups = remove_current_folds_from_equal_rule_groups(current_folds, equal_rule_groups)
#     # print(f"{future_rule2rns=}")
#     # equal_rule_dict = [rule_sizes(i) for i in equal_rule_groups]

#     # print("\tEliminated folds/rules that enable subsets of others' reactions.")
#     # print("\t\t%i folds available for the NEXT ITERATION"%len(set([f for g in equal_rule_groups for r in g for f in r])))
#     # print("\t\t%i rules available for the NEXT ITERATION"%len(equal_rule_groups))
#     return [rule_sizes(i) for i in equal_rule_groups]

def next_iter_possible_rules(current_fold2s, scope_rule2rn, remaining_rules, current_rns):
    """
    Returns a list of equal rule group dictionaries, keyed by rule size.

    :param current_folds: collection of current folds
    :param scope_rule2rn: dict of rule:rns (should be for scope)
    :return: a list of equal rule group dictionaries, keyed by rule size
            [
                {1:[rule1,rule2,rule2],3:[rule4,rule5]}, #dict of equal rules, 3 with len=1, 2 with len=3
                {1:[rule7,rule8],3:[rule9]}              
            ]
    """

    ## Need to run these two calls every iteration of the fold expansion
    # future_rule2rns = rule2nextrns(current_folds, scope_rule2rn) ## This is just finnding rules from reamaining folds with current reactions appended to all values
    future_rule2rns = {k:(v | current_rns) for k,v in remaining_rules.items()}
    equal_rule_groups = create_equal_rule_groups(future_rule2rns) ## discards any rules which are strict subsets of others...
    # equal_rule_groups = sort_equal_rule_groups(equal_rule_groups)
    # new_equal_rule_groups = []
    # for superset,subset in equal_rule_groups.items():
    #     i = frozenset([rule-current_folds for rule in superset])
    #     j = frozenset([rule-current_folds for rule in subset])
    #     new_equal_rule_groups.append({i:j})
    ## Remove current folds from each equivilent rule group
    for i in equal_rule_groups:
        i.equal_supersets = [frozenset(rule - current_folds) for rule in i.equal_supersets]
        i.equal_subsets = [frozenset(rule - current_folds) for rule in i.equal_subsets]

        i.equal_supersets_by_size = rule_sizes(i.equal_supersets)
        i.equal_subsets_by_size = rule_sizes(i.equal_subsets)
    
    equal_rule_groups = set(list(equal_rule_groups)) ## gets rid of duplicates
    # print(f"{future_rule2rns=}")
    # equal_rule_dict = [rule_sizes(i) for i in equal_rule_groups]

    # print("\tEliminated folds/rules that enable subsets of others' reactions.")
    # print("\t\t%i folds available for the NEXT ITERATION"%len(set([f for g in equal_rule_groups for r in g for f in r])))
    # print("\t\t%i rules available for the NEXT ITERATION"%len(equal_rule_groups))
    # return [(rule_sizes(i), rulesizes(j)) for i,j in equal_rule_groups.items()]
    return equal_rule_groups

def maxreactions(r_effects):
    """
    Returns the rule enabling the most reactions

    :param r_effects: a dict of rule:{rule2rn:... , cpds:... , rns:...} 
                      denoting the outcome of adding each rule.
    :return: rule enabling the most reactions
    """
    k_vcount = {k:len(v["rns"]) for k,v in r_effects.items()}
    k_vcount = dict(sorted(k_vcount.items())) ## Sort for reproduceability
    return max(k_vcount, key = k_vcount.get)

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

# def free_rules(current_folds, scope_rules2rn):
#     """
#     Returns rules that weren't explicity added, yet whose reactions are already enabled.
    
#     This can occur for example if a rule is a subset of another rule which was selected.

#     Not all folds in free_rules are free
#     e.g. F1,F2 -> R1, F3 -> R1, F1 -> R7
#     If F3 already discovered, even though the first rule is redundent, it can't be "free"
#     because F1 additionally enables a new reaction (R7) which hasn't yet been discovered.

#     :param current_folds: collection of current folds
#     :param scope_rules2rn: dict of scope rules2rn
#     :return: a set of rules whose folds are not part of current_folds, yet whose reactions
#                 are all already enabled.
#     """
#     current_rule2rn = subset_rule2rn_from_folds(current_folds, scope_rules2rn)
#     current_rns = set([rn for v in current_rule2rn.values() for rn in v])
#     return {k for k,v in scope_rules2rn.items() if (v <= current_rns) and not (k <= current_folds)}


########################################################################################################################
########################################################################################################################
class EquivilentRule:
    """
    """
    ## if the equivilentrule.supersets are the same don't add the second one to the set of all equivilent rules

    def __init__(self, equal_supersets, equal_subsets):
        """
        """

        self.equal_supersets = frozenset(equal_supersets)
        self.equal_subsets = frozenset(equal_subsets)
        self.equal_supersets_by_size = None
        self.equal_subsets_by_size = None

    def __eq__(self, other):
        return self.equal_supersets == other.equal_supersets

    def __hash__(self):
        return hash(equal_supersets)

class FoldRules:
    """
    Stores how fold rules and reactions map to each other, and which reactions are independent of folds.

    Invariant to seed compounds, seed folds, etc. This is the universe of possible fold rules. 
    """

    def __init__(self, rn2rules, fold_independent_rns):
        """
        Construct fold rules from an `rn2rules` dictionary, and a set of `fold_independent_rns`.

        :param rn2rules: dict of rn:{rules}, where rules are frozensets of folds
        :param fold_independent_rns: a set of reactions which can occur without folds
        """

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rule2rns = rule2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rule2rns.keys() for i in fs]) ## all folds
        self.fold_independent_rns = fold_independent_rns

class FoldMetabolism:
    ## Rename to something like FoldScope?
    """
    A class to do fold expansion from a metabolism, foldrules, and seeds.
    """

    def __init__(self, metabolism, foldrules, preexpansion=False):
        """
        Calculates expansion scope after seed compounds are defined. 

        :param metabolism: a GlobalMetabolicNetwork object that defines the compound/reaction rules of network expansion
        :param foldrules: a FoldRules object that defines the fold rules of fold expansion
        :kwarg preexpansion: NOT YET IMPLEMENTED -- but the idea is to allow n steps of regular network expansion before beginning fold expansion
        """
        
        self._m = metabolism ## GlobalMetabolicNetwork object
        self._f = foldrules # FoldRules object

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
            self.scope_rn2rules = {k:v for k,v in self.f.rn2rules.items() if k in self.scope_rns}
            self.scope_rules2rn = rule2rn(self.scope_rn2rules)
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

    # def loop_through_rules(self, current_folds, current_cpds, current_rns):
    #     """
    #     Loops through all remaining rules in the scope and returns a dict 
    #         of rule:{rule2rn:... , cpds:... , rns:...} denoting the
    #         outcome of adding each rule.

    #     :param current_folds: collection of current folds
    #     :param current_cpds: collection of current cpds
    #     :param current_rns: collection of current rns
    #     :return: - a dict of rule:{rule2rn:... , cpds:... , rns:...} 
    #                denoting the outcome of adding each rule.
    #              - n_rules_checked for metadata purposes
    #              - n_equal_rule_groups for metadata purposes
    #     """

    #     equal_rule_dict = next_iter_possible_rules(current_folds2, self.scope_rules2rn)
    #     if len(equal_rule_dict) == 0:
    #         print("len(equal_rule_dict) == 0, reffects will be empty")
    #         r_effects = dict()

    #     rule_sizes = set(list([i for d in equal_rule_dict for i in d.keys()]))
    #     n_rules_checked = 0 # for metadata
    #     for rsize in sorted(rule_sizes):
    #         r_effects = dict()
    #         for d in equal_rule_dict:
    #             if rsize in d:
    #                 # only look at first rule among equals
    #                 rule = d[rsize][0]
    #                 _fdict = dict()
    #                 _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(rule, current_folds, current_cpds)
    #                 r_effects[rule] = _fdict
    #                 n_rules_checked+=1

    #         ## Don't look for longer rules if shorter rules enable new reactions
    #         if len(r_effects) > 0:
    #             break

    #     ## r_effects now only fills will rules that actually end up adding reactions
    #     return r_effects, n_rules_checked, len(equal_rule_dict)

    def loop_through_rules(self, current_folds, current_cpds, current_rns):
        """
        Loops through all remaining rules in the scope and returns a dict 
            of rule:{rule2rn:... , cpds:... , rns:...} denoting the
            outcome of adding each rule.

        :param current_folds: collection of current folds
        :param current_cpds: collection of current cpds
        :param current_rns: collection of current rns
        :return: - a dict of rule:{rule2rn:... , cpds:... , rns:...} 
                   denoting the outcome of adding each rule.
                 - n_rules_checked for metadata purposes
                 - n_equal_rule_groups for metadata purposes
        """

        equal_rule_dict = next_iter_possible_rules(current_folds, self.scope_rules2rn)
        # if len(equal_rule_dict) == 0:
        #     print("len(equal_rule_dict) == 0, reffects will be empty")
        #     r_effects = dict()

        ## LEFT OFF HERE ON NOV 25--NEED TO LOOP THROUGH MY EQUAL_RULE_DICT OF OBJECTS TO REWRITE THE NEXT BIT HERE
        ##      TAKING INTO ACCOUNT THAT I NEED TO CHANGE HOW I BREAK SINCE ALL R_EFFECTS WILL BE >0 NOW, AND
        ##      THAT I HAVE TO SEARCH THE EQUAL_RULE_DICT_SUBSETS FOR WHATEVER SUPERSETS ARE FOUND TO GENERATE THE MAX
        ##      NUMBER OF REACTIONS. Only have to search subsets which have the same rulesize as the discovered superset though.
        ## It could also be that some subsets have shorter rule lengths than the supersets, which complicates things. 

        superset_rule_sizes = set(list([i for d in equal_rule_dict for i in d[0].keys()]))
        n_rules_checked = 0 # for metadata
        for rsize in sorted(superset_rule_sizes):
            r_effects = dict()
            for d_superset, d_subset in equal_rule_dict:
                if rsize in d_superset:
                    rule = random.choice(d_superset[rsize][0]) ## chose a random rule from equivilent rules
                    _fdict = dict()
                    _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(rule, current_folds, current_cpds)
                    r_effects[rule] = _fdict
                    n_rules_checked+=1

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(r_effects) > 0:
                break

        ## r_effects now only fills will rules that actually end up adding reactions
        return r_effects, n_rules_checked, len(equal_rule_dict)

    def select_next_rule_or_fold(self, current_folds, current_cpds, current_rns, remaining_folds, remaining_rules, algorithm="maxreactions"):#rselect_func=maxreactions):
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
            return next_fold, _fdict, 1, None

        elif algorithm == "randomrule":
            rule_size_dict = rule_sizes(remaining_rules)
            next_rule = random.choice(rule_size_dict[min(rule_size_dict)])
            _fdict = dict()
            _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(next_fold, current_folds, current_cpds)
            return next_rule, _fdict, 1, None
            
        elif algorithm == "maxreactionsupersets":
            r_effects, n_rules_checked, n_equal_rule_groups = self.loop_through_rules(current_folds, current_cpds, current_rns)
            # if len(r_effects) == 0:
            #     next_rule = None
            #     r_effects[next_rule] = {"cpds":deepcopy(current_cpds), "rns":deepcopy(current_rns)}
            # else:
            #     next_rule = rselect_func(r_effects)
            next_rule = rselect_func(r_effects)
            return next_rule, r_effects[next_rule], n_rules_checked, n_equal_rule_groups

        elif algorithm == "maxreactions":
            r_effects, n_rules_checked, n_equal_rule_groups = self.loop_through_rules(current_folds, current_cpds, current_rns)
    
    def rule_order(self, algorithm="randomfold"):
        """
        Determine the ordering of all rules/folds.

        :kwarg track_free_rules: if True, add rules/folds to the iteration_dict that
                           weren't selected, but whose reactions are all
                           already enabled.

                           (Probably I want to just remove this kwarg and
                            always have this enabled...)
        """
        valid_algorithms = ["randomfold", "randomrule", "maxreactions"]
        if algorithm.lower() not in valid_algorithms:
            raise ValueError("algorithm must be one of %s"%valid_algorithms)

        if (self.seed_cpds == None) or (self.seed_folds == None):
            raise ValueError("self.seed_cpds and self.seed_folds must not be None")
        
        ## Data to store meta information about run
        metadict = {
            "runtime": dict(),
            "freefolds": dict(),
            "n_rulegroups_in_iteration":dict(),
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
        metadict["n_rulegroups_in_iteration"][iteration] = 0
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
            next_rule, fdata, n_rules_checked, n_equal_rule_groups = self.select_next_rule_or_fold(current["folds"], current["cpds"], current["rns"], remaining_folds, remaining_rules, algorithm)
            remaining_folds = (remaining_folds - set(next_rule))
            remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if len(k & remaining_folds)>0}
            # remaining_rules = {k:v for k,v in self.scope_rules2rn.items() if k not in subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn)}

            ## Stop conditions
            if len(remaining_folds) == 0:
                keepgoing = False
            if next_rule == None:
                keepgoing == False
            elif (fdata["cpds"] == current["cpds"]) and (fdata["rns"] == current["rns"]):
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
            metadict["n_rulegroups_in_iteration"][iteration] = n_equal_rule_groups
            metadict["n_rules_checked"][iteration] = n_rules_checked
            metadict["max_n_remaining_folds"][iteration] = len(remaining_folds)
            metadict["max_n_remaining_rules"][iteration] = len(self.scope_rules2rn) - len(subset_rule2rn_from_folds(current["folds"], self.scope_rules2rn))

        return current, iteration_dict, metadict

def example_main():
    ## Load metabolism
    metabolism_path = PurePath(asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_path)
    
    ## Load fold rules
    rn2rules_db_path = PurePath("data", "rn2fold", "_db.pkl")
    rn2rules_db = pd.read_pickle(rn2rules_db_path)

    rn2rules_path = rn2rules_db.iloc[4]["OUTPUT_PATH"]
    rn2rules = pd.read_pickle(rn2rules_path)

    fold_independent_rns = set()
    foldrules = nf.FoldRules(rn2rules, fold_independent_rns)

    ## Inititalize fold metabolism
    fm = nf.FoldMetabolism(metabolism, foldrules)

    seed_cpds_path = PurePath("data", "josh", "seed_set.csv")
    fm.seed_cpds = set((pd.read_csv(seed_cpds_path)["ID"]))
    fm.seed_folds = set(['spontaneous'])

    ## Run fold expansion
    current, iteration_dict, metadict = fm.rule_order()
    