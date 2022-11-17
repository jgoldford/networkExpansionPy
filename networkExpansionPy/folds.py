import pandas as pd
from pathlib import PurePath, Path
from copy import deepcopy
import timeit
from pprint import pprint
from collections import Counter

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

def subset_rule2rn(folds, rule2rn):
    """
    Returns a dictionary of rules:rns enabled by folds
    
    :param folds: collection of folds used to subset the rules dict
    :param rules: dict of rule:rns mappings to subset from
    :return: dictionary of rules:rns enabled by folds
    """
    return {k:v for k,v in rule2rn.items() if k <= set(folds)}

def rule2nextrns(current_folds, rule2rn):
    """
    Returns a dictionary of rules:rns only for rules enabling reactions that 
        are undiscovered by current_folds. Each rule will map to all reactions
        enabled after that rule is added.

    :param current_folds: collection of folds to compare to
    :param rule2rn: dict of rule:rns mappings to subset from, and then compare to
    """
    current_rule2rn = subset_rule2rn(current_folds, rule2rn)
    current_rns = set([rn for v in current_rule2rn.values() for rn in v])
    return {k:(v | current_rns) for k,v in rule2rn.items() if not v <= current_rns}

def create_equal_rule_groups(rule2rn):
    """
    Returns a set of equivilent rule collections.

    Collections and the rules they contain are each frozensets.

    :param rule2rn: dict of rule:rns to create equal rule groups from
    :return: set of frozensets (i.e. groups of equivilent rules) of frozensets (i.e. rules) 
    """

    strictsubset_ks = set()
    equal_groups = set()
    for k, v in rule2rn.items():
        equal_ks = {k}
        for k2, v2 in rule2rn.items():
            if k != k2:
                if v2 <= v:
                    if v2!=v:
                        strictsubset_ks.add(k2)
                    else: # v2==v
                        equal_ks.add(k2)
        
        equal_groups.add(frozenset(equal_ks))

    ## exclude groups which include any strict subset rules
    return {i for i in equal_groups if not i & strictsubset_ks}

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

def remove_current_folds_from_equal_rule_groups(current_folds, rule_groups):
    """ 
    Returns a list of rule groups with current folds removed.

    Preserves ordering from `rule_groups`.

    :param current_folds: collection of current folds
    :param rule_groups: a collection of equivilent rule groups
    :return: a list of rule groups with current folds removed
    """

    new_rule_groups = []
    for rg in rule_groups:
        new_group = []
        # this construction allows us to preserve the ordering from create_equal_rule_groups
        for i in rg:
            new_rule = frozenset(i) - current_folds
            if len(new_rule) > 0:
                new_group.append(new_rule)
        if len(new_group) > 0:
            new_rule_groups.append(new_group)

    return new_rule_groups

def next_iter_possible_rules(current_folds, rule2rn):
    """
    Returns a list of equal rule group dictionaries, keyed by rule size.

    :param current_folds: collection of current folds
    :param rule2rn: dict of rule:rns (should be for scope)
    :return: a list of equal rule group dictionaries, keyed by rule size
            [
                {1:[rule1,rule2,rule2],3:[rule4,rule5]}, #dict of equal rules, 3 with len=1, 2 with len=3
                {1:[rule7,rule8],3:[rule9]}              
            ]
    """

    ## Need to run these two calls every iteration of the fold expansion
    future_rule2rns = rule2nextrns(current_folds, rule2rn)
    equal_rule_groups = create_equal_rule_groups(future_rule2rns)
    equal_rule_groups = sort_equal_rule_groups(equal_rule_groups)
    equal_rule_groups = remove_current_folds_from_equal_rule_groups(current_folds, equal_rule_groups)
    # print(f"{future_rule2rns=}")
    # equal_rule_dict = [rule_sizes(i) for i in equal_rule_groups]

    # print("\tEliminated folds/rules that enable subsets of others' reactions.")
    # print("\t\t%i folds available for the NEXT ITERATION"%len(set([f for g in equal_rule_groups for r in g for f in r])))
    # print("\t\t%i rules available for the NEXT ITERATION"%len(equal_rule_groups))
    return [rule_sizes(i) for i in equal_rule_groups]

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

########################################################################################################################
########################################################################################################################
class FoldRules:

    def __init__(self, rn2rules, fold_independent_rns):

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rule2rns = rule2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rule2rns.keys() for i in fs]) ## all folds
        self.fold_independent_rns = fold_independent_rns

class FoldMetabolism:
    """
    A re-write of FoldMetabolism to handle doing network expansion with folds
    """

    def __init__(self, metabolism, foldrules, preexpansion=False):
        # load the data
        self._m = metabolism ## Metabolism object
        self._f = foldrules

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
            print("Seed compounds updated. Calculating seed's scope (this is only done once) ... ")
            self._seed_cpds = seed_cpds 
            self.scope_cpds, self.scope_rns = self.calculate_scope(self._seed_cpds)
            self.scope_rn2rules = {k:v for k,v in self.f.rn2rules.items() if k in self.scope_rns}
            self.scope_rules2rn = rule2rn(self.scope_rn2rules)
            self.scope_folds = set([i for fs in self.scope_rules2rn.keys() for i in fs])
            print("... done.")
            # print("\tEliminated folds/rules outside of scope or fold independent reactions.") ## scope only includes reactions in folds, and reacitons explicitly input as independent of folds
            # print("\t\t%i folds available in RUN"%len(self.scope_folds))
            # print("\t\t%i rules available in RUN"%len(self.scope_rules2rn))

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

    def effect_per_rule(self, rule, current_folds, current_cpds):
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
        potential_rule2rns = subset_rule2rn(potential_fold_set, self.scope_rules2rn)
        cx,rx = self.fold_expand(potential_rule2rns, current_cpds)

        return potential_rule2rns, cx, rx

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

        rule_sizes = set(list([i for d in equal_rule_dict for i in d.keys()]))
        # print(f"{rule_sizes=}")
        # print(f"{len(equal_rule_dict)=}")
        # print("lens of rules in equal_rule_dict: ", Counter([k for d in equal_rule_dict for k,v in d.items()]))
        # print("n rules remaining: ", len(equal_rule_dict))
        # print("")
        # print(f"{equal_rule_dict=}")
        n_rules_checked = 0 # for metadata
        for rsize in sorted(rule_sizes):
            r_effects = dict()
            for d in equal_rule_dict:
                if rsize in d:
                    # only look at first rule among equals
                    rule = d[rsize][0]
                    _fdict = dict()
                    _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule(rule, current_folds, current_cpds)
                    if len(_fdict["rns"] - current_rns) > 0:
                        r_effects[rule] = _fdict
                    n_rules_checked+=1

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(r_effects) > 0:
                break

        ## r_effects now only fills will rules that actually end up adding reactions
        return r_effects, n_rules_checked, len(equal_rule_dict)

    def select_next_rule(self, current_folds, current_cpds, current_rns, rselect_func=maxreactions):
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
        r_effects, n_rules_checked, n_equal_rule_groups = self.loop_through_rules(current_folds, current_cpds, current_rns)
        next_rule = rselect_func(r_effects)
        return next_rule, r_effects[next_rule], n_rules_checked, n_equal_rule_groups
    
    def free_rules(self, current_rns, current_folds):
        """
        Returns rules that weren't explicity added, yet whose reactions are already enabled.
        
        This can occur for example if a rule is a subset of another rule which was selected.

        :current_rns: collection of current reactions
        :current_folds: collection of current folds
        :return: a set of rules whose folds are not part of current_folds, yet whose reactions
                 are all already enabled.
        """
        return {k for k,v in self.scope_rules2rn.items() if (v <= current_rns) and not (k <= current_folds)}

    def rule_order(self, free_rules=True):
        """
        Determine the ordering of all rules/folds.

        :kwarg free_rules: if True, add rules/folds to the iteration_dict that
                           weren't selected, but whose reactions are all
                           already enabled.

                           (Probably I want to just remove this kwarg and
                            always have this enabled...)
        """

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
        
        ################
        ## ITERATION 0
        ## Avoid updating folds on the 0th iteration since they don't apply until iteration=1
        iteration_dict = update_iteration_dict(iteration_dict, {k:v for k,v in current.items() if k!="folds"}, iteration)
        iteration+=1
        start = timeit.default_timer()

        # print("iteration dict INIT: ", iteration_dict)
        # print("current dict INIT: ", current)

        ################
        ## ITERATION 1 (using only seed folds and fold independent reactions)
        init_rules2rn = subset_rule2rn(current["folds"], self.scope_rules2rn)
        current["cpds"], current["rns"] = self.fold_expand(init_rules2rn, current["cpds"])
        ## Add free folds to current dict
        free_folds = {i for fs in self.free_rules(current["rns"], current["folds"]) for i in fs}
        # print("current folds: ", current["folds"])
        # print("free folds: ", free_folds)
        if free_rules == True: ## Append free_folds to data dict
            current["folds"] = (current["folds"] | free_folds)
        remaining_folds = (self.scope_folds - current["folds"] - free_folds) ## Remove the free folds from the remaining folds regardless
        iteration_dict = update_iteration_dict(iteration_dict, current, iteration)

        ## Update metadata
        metadict["runtime"][iteration] = timeit.default_timer() - start
        metadict["freefolds"][iteration] = free_folds
        metadict["n_rulegroups_in_iteration"][iteration] = 0
        metadict["n_rules_checked"][iteration] = 0
        metadict["max_n_remaining_folds"][iteration] = len(remaining_folds)
        metadict["max_n_remaining_rules"][iteration] = len(self.scope_rules2rn) - len(subset_rule2rn(current["folds"], self.scope_rules2rn))
        # print("\nITERATION: ", iteration)
        # for k, d in metadict.items():
        #     if k!="freefolds":
        #         print(k,": ",d[iteration])

        ## Needed in case expansion not possible at all
        if len(remaining_folds) > 0:
            keepgoing = True
        else:
            keepgoing = False

        # print("iteration dict after first expansion: ", iteration_dict)
        # print("current dict after first expansion: ", current)
        ################
        ## ITERATION 2+
        while keepgoing:
            start = timeit.default_timer()
            iteration += 1
            # print("\nITERATION: ", iteration)
            # print("maximum n folds remaining: ", len(remaining_folds))
            # print("maximum remaining_folds:\n", remaining_folds)
            # print("subset_rule2rn: ", subset_rule2rn(remaining_folds, self.scope_rules2rn))
            next_rule, fdata, n_rules_checked, n_equal_rule_groups = self.select_next_rule(current["folds"], current["cpds"], current["rns"])
            # print("next fold: ", next_rule)
            # print("fdata:")
            # pprint(fdata)
            free_folds = {i for fs in self.free_rules(fdata["rns"], (current["folds"] | set(next_rule))) for i in fs}
            # print("free folds: ", free_folds)
            remaining_folds = (remaining_folds - set(next_rule) - free_folds)

            if len(remaining_folds) == 0:
                keepgoing = False

            ## Stop conditions
            if (fdata["cpds"] == current["cpds"]) and (fdata["rns"] == current["rns"]):
                keepgoing = False    

            else:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                if free_rules == True: ## Append free_folds to data dict
                    current["folds"] = (current["folds"] | set(next_rule) | free_folds)
                else:
                    current["folds"] = (current["folds"] | set(next_rule))
                current["cpds"] = fdata["cpds"]
                current["rns"] = fdata["rns"]
                
                ## Store when cpds and rns appear in the expansion
                iteration_dict = update_iteration_dict(iteration_dict, current, iteration)

            ## Update metadata
            metadict["runtime"][iteration] = timeit.default_timer() - start
            metadict["freefolds"][iteration] = free_folds
            metadict["n_rulegroups_in_iteration"][iteration] = n_equal_rule_groups
            metadict["n_rules_checked"][iteration] = n_rules_checked
            metadict["max_n_remaining_folds"][iteration] = len(remaining_folds)
            metadict["max_n_remaining_rules"][iteration] = len(self.scope_rules2rn) - len(subset_rule2rn(current["folds"], self.scope_rules2rn))
            # print("\nITERATION: ", iteration)
            # for k, d in metadict.items():
            #     if k!="freefolds":
            #         print(k,": ",d[iteration])

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
    