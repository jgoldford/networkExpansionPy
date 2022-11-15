from operator import index
import numpy as np
import pandas as pd
from pathlib import PurePath, Path
from copy import copy, deepcopy
import timeit
from pprint import pprint
from collections import Counter

asset_path = PurePath(__file__).parent / "assets"

# define a function that determins if a reaction rule (x) is feasible with foldSet
def singlerule2rn(foldSet,x):
    if x.issubset(foldSet):
    #if foldSet.issubset(x)
        return True
    else:
        return False

# define a function that returns a list of reactions that are feasible with foldSet.  rules_sub contains all the fold set rules for all reactions
def dffolds2rn(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: singlerule2rn(foldSet,x))
    y = rules_sub[feasible]
    if len(y) > 0:
        rns = rules_sub[feasible].rn.unique().tolist()
    else:
        rns = []
    return rns

# same as above but return only the set of rules that are feasible
def dffolds2rules(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: singlerule2rn(foldSet,x))
    rule_df = rules_sub[feasible]
    return rule_df


# define a fold rules object, which simply contains fold sets that enable reactions
class FoldRules:
    
    def __init__(self, path = PurePath("ecode","ecod2rn.ec3.07Feb2021.csv")):
        # load the data
        # self.allrules = setRules(path=path) ## This should be all rules, every step of expansion will subset
        self.rules = None; ## This should be current rules
        self.foldrns = None ## reactions enabled by all fold rules
        self.folds = None; ## all folds

        ## probably my convention of how i store these values will change if I read in a dataframe vs my pickled dictionary
        ## - current rules should also be part of the fold expansion object
        ## - probably fold rules don't even need to be an object themselves


    def copy(self):
        return deepcopy(self)

    def setRules(self,path = PurePath("ecode","ecod2rn.ec3.07Feb2021.csv")):
        rules = pd.read_csv(asset_path / path, index_col=[0]).reset_index(drop=True)
        
        self.foldrns = rules.rn.unique().tolist()
        folds =  rules['rule'].apply(lambda x: set(x.split('_')))
        folds = [item for sublist in folds for item in sublist]
        self.folds = list(set(folds))
        rules['fold_sets']= rules['rule'].apply(lambda x: set(x.split('_')))
        self.rules = rules
        
    def folds2reactions(self,foldSet):
        # make sure foldSet is of "set" type
        if type(foldSet) != set:
            if type(foldSet) != list:
                foldSet = set([foldSet])
            else:
                foldSet = set(foldSet)
                
        rns = dffolds2rn(self.rules,foldSet)
        return rns
    
    def removeFolds(self,folds_remove):
        if type(folds_remove) != set:
            if type(folds_remove) != list:
                folds_remove = set([folds_remove])
            else:
                folds_remove = set(folds_remove)
                
        rulesRemoval = self.rules.fold_sets.apply(lambda x: folds_remove.issubset(x)) 
        self.rules = self.rules[~rulesRemoval]
        folds_remove_list = list(folds_remove)
        self.folds = [x for x in self.folds if x not in folds_remove_list]

#function that peforms fold expansion
def fold_expansion(metabolism,foldRules,fold_set,cpds_set,rxns_seed):
    rn_list = foldRules.folds2reactions(fold_set)
    rn_list = metabolism.rxns2tuple(rn_list) + list(rxns_seed)
    cx,rx = metabolism.expand(cpds_set,reaction_mask=rn_list)
    return cx,rx,rn_list

########################################################################################################################
########################################################################################################################
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

def rule2rn_enabling_new_rn(current_folds, rule2rn):
    """
    Returns a dictionary of rules:rns only for rules enabling reactions that 
        are undiscovered by current_folds.

    :param current_folds: collection of folds to compare to
    :param rule2rn: dict of rule:rns mappings to subset from, and then compare to
    """
    current_rule2rn = subset_rule2rn(current_folds, rule2rn)
    current_rns = set([rn for v in current_rule2rn.values() for rn in v])
    return {k:(v | current_rns) for k,v in rule2rn.items() if not v <= current_rns}

def create_equal_rule_groups(rule2rn):
    """
    Returns a sorted list of equivilent rule collections.

    Collections and list are each sorted.

    :param rule2rn: dict of rule:rns to create equal rule groups from
    :return: a sorted list of equivilent rule collections (each sorted)
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
    equal_groups = {i for i in equal_groups if not i & strictsubset_ks}

    ## transform to lists of lists of rules for the purposes of sortability/reproducability
    return [sorted([sorted(i) for i in group]) for group in equal_groups]
########################################################################################################################

class GlobalFoldNetwork:

    def __init__(self, rn2rules, fold_independent_rns):

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rule2rns = rule2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rule2rns.keys() for i in fs]) ## all folds
        print("GlobalFoldNetwork initialized\n%i folds available in RUN"%(len(self.folds)))
        self.fold_independent_rns = fold_independent_rns

class FoldMetabolism:
    """
    A re-write of FoldMetabolism to handle doing network expansion with folds
    """

    def __init__(self, metabolism, foldnet, preexpansion=False):
        # load the data
        self._m = metabolism ## Metabolism object
        self._f = foldnet

        self.seed_folds = None
        self._seed_cpds = None

        self.scope_cpds = None
        self.scope_rns = None
        self.scope_folds = None
        self.scope_rules2rn = None
        self.scope_rn2rules = None

    ## Disallow changing metabolism or foldnet after initialization b/c no setter
    @property 
    def metabolism(self):
        return self._m
    
    @property 
    def foldnet(self):
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
            self.scope_rn2rules = {k:v for k,v in self._f.rn2rules.items() if k in self.scope_rns}
            self.scope_rules2rn = rule2rn(self.scope_rn2rules)
            self.scope_folds = set([i for fs in self.scope_rules2rn.keys() for i in fs])
            print("...done.")
            print("Folds in RUN reduced to overlap with scope from fold reactions and fold independent reactions\n%i folds available in RUN"%len(self.scope_folds)) ## scope only includes reactions in folds, and reacitons explicitly input as independent of folds

        else:
            pass 
                
    def calculate_scope(self, seed_cpds):
        """
        Calculate the scope of the seed_cpds with all reactions enabled by the global fold network (including fold_independent reactions)

        :param seed_cpds: collection of compound ids
        :return: set of compounds and set of reactions in the scope
        """
        rn_tup_set = set(self._m.rxns2tuple((self._f.rns | self._f.fold_independent_rns)))
        scope_cpds, scope_rns = self._m.expand(seed_cpds, reaction_mask=rn_tup_set)
        return set(scope_cpds), set([i[0] for i in scope_rns])

    def organize_set_elements_by_size(self, myset):
        element_lengths = dict()
        for i in myset:
            if not len(i) in element_lengths:
                element_lengths[len(i)] = [i]
            else:
                element_lengths[len(i)].append(i)
        return element_lengths

    def organize_equal_rule_groups_by_size(self, listoflists):
        ## each inner list should turn into a dict organized by size
        outer_list = []
        for i in listoflists:
            outer_list.append(self.organize_set_elements_by_size(i))
        return outer_list

    def remove_current_folds_from_equal_rule_groups(self, current_folds, rule_groups):

        new_rule_groups = []
        for rg in rule_groups:
            new_group = []
            # this construction allows us to preserve the ordering from create_equal_rule_groups
            for i in rg:
                new_rule = frozenset(i) - current_folds
                if len(i) > 0:
                    new_group.append(new_rule)
            if len(new_group) > 0:
                new_rule_groups.append(new_group)

        return new_rule_groups
    
    def next_iter_possible_rules(self, current_folds):

        ## Need to run these two calls every iteration of the fold expansion
        future_rule2rns = rule2rn_enabling_new_rn(current_folds, self.scope_rules2rn)
        print(f"{future_rule2rns=}")
        equal_rule_groups_that_arent_subsets = create_equal_rule_groups(future_rule2rns)
        equal_rule_dict = self.organize_equal_rule_groups_by_size(self.remove_current_folds_from_equal_rule_groups(current_folds, equal_rule_groups_that_arent_subsets))

        # print("-> Folds whose rules correspond to reactions which are subsets of one another in the next NEXT ITERATION removed\n-> ... %i folds available for the NEXT ITERATION"%len(filtered_folds_to_expand))
        return equal_rule_dict #filtered_folds_to_expand

    def fold_expand(self, metabolism, rule2rns, fold_independent_rns, cpds):
        """Doesn't use self"""
        fold_rns = set([i for s in rule2rns.values() for i in s])
        rn_tup_set = set(metabolism.rxns2tuple(fold_rns | fold_independent_rns))
        cx,rx = metabolism.expand(cpds,reaction_mask=rn_tup_set)
        return set(cx), set([i[0] for i in rx])

    def effect_per_rule(self, rule, current_folds, current_cpds):
        """
        self.scope_rules2rn
        self._m
        """

        potential_fold_set = (current_folds | set(rule))
        potential_rules2rn = subset_rule2rn(potential_fold_set, self.scope_rules2rn)
        cx,rx = self.fold_expand(self._m, potential_rules2rn, self._f.fold_independent_rns, current_cpds)

        return potential_rules2rn, cx, rx #set(cx), set(rx)

    def loop_through_rules(self, current_folds, current_cpds, current_rns):

        equal_rule_dict = self.next_iter_possible_rules(current_folds)

        rule_sizes = set(list([i for d in equal_rule_dict for i in d.keys()]))
        print(f"{rule_sizes=}")
        print(f"{len(equal_rule_dict)=}")
        print("lens of rules in equal_rule_dict: ", Counter([k for d in equal_rule_dict for k,v in d.items()]))
        print("n rules remaining: ", len(equal_rule_dict))
        print("")
        print(f"{equal_rule_dict=}")

        for fsize in sorted(rule_sizes):
            f_effects = dict()
            for d in equal_rule_dict:
                if fsize in d:
                    # only look at first rule among equals
                    rule = d[fsize][0]
                    _fdict = dict()
                    _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule(rule, current_folds, current_cpds)
                    if len(_fdict["rns"] - current_rns) > 0:
                        f_effects[rule] = _fdict

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(f_effects) > 0:
                break

        ## f_effects now only fills will rules that actually end up adding reactions
        return f_effects

    def maxreactions(f_effects):
        """Doesn't use self"""
        k_vcount = {k:len(v["rns"]) for k,v in f_effects.items()}
        k_vcount = dict(sorted(k_vcount.items())) ## Sort for reproduceability
        return max(k_vcount, key = k_vcount.get)

    def select_next_rule(self, current_folds, current_cpds, current_rns, fselect_func=maxreactions):
        """Doesn't use self"""
        f_effects = self.loop_through_rules(current_folds, current_cpds, current_rns)
        next_rule = fselect_func(f_effects)
        return next_rule, f_effects[next_rule]

    def update_iteration_dict(self, iteration_dict, current, iteration):
        for dtype, ids in current.items():
            for i in ids:
                if i not in iteration_dict[dtype]:
                    iteration_dict[dtype][i] = iteration
        return iteration_dict
    
    def free_rules(self, current_rns, current_folds):
        return {k for k,v in self.scope_rules2rn.items() if (v <= current_rns) and not (k <= current_folds)}

    def rule_order(self, free_rules=True):

        if (self.seed_cpds == None) or (self.seed_folds == None):
            raise ValueError("self.seed_cpds and self.seed_folds must not be None")
        
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
        
        ## Avoid updating folds on the 0th iteration since they don't apply until iteration=1
        iteration_dict = self.update_iteration_dict(iteration_dict, {k:v for k,v in current.items() if k!="folds"}, iteration)
        iteration+=1

        # print("iteration dict INIT: ", iteration_dict)
        # print("current dict INIT: ", current)

        ## First expansion (using only seed folds and fold independent reactions)
        init_rules2rn = subset_rule2rn(current["folds"], self.scope_rules2rn)
        current["cpds"], current["rns"] = self.fold_expand(self._m, init_rules2rn, self._f.fold_independent_rns, current["cpds"])
        ## Add free folds to current dict
        free_folds = {i for fs in self.free_rules(current["rns"], current["folds"]) for i in fs}
        print("current folds: ", current["folds"])
        print("free folds: ", free_folds)
        if free_rules == True: ## Append free_folds to data dict
            current["folds"] = (current["folds"] | free_folds)
        remaining_folds = (self.scope_folds - current["folds"] - free_folds) ## Remove the free folds from the remaining folds regardless
        iteration_dict = self.update_iteration_dict(iteration_dict, current, iteration)

        ## Needed in case expansion not possible at all
        if len(remaining_folds) > 0:
            keepgoing = True
        else:
            keepgoing = False

        # print("iteration dict after first expansion: ", iteration_dict)
        # print("current dict after first expansion: ", current)

        while keepgoing:
            start = timeit.default_timer()
            iteration += 1
            print("\nITERATION: ", iteration)
            print("maximum n folds remaining: ", len(remaining_folds))
            print("maximum remaining_folds:\n", remaining_folds)
            print("subset_rule2rn: ", subset_rule2rn(remaining_folds, self.scope_rules2rn))
            next_rule, fdata = self.select_next_rule(current["folds"], current["cpds"], current["rns"])
            print("next fold: ", next_rule)
            print("fdata:")
            pprint(fdata)
            exec_time = timeit.default_timer() - start
            print("iteration runtime: ", exec_time)
            free_folds = {i for fs in self.free_rules(fdata["rns"], (current["folds"] | set(next_rule))) for i in fs}
            print("free folds: ", free_folds)
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
                iteration_dict = self.update_iteration_dict(iteration_dict, current, iteration)

        return current, iteration_dict

def example_main():
    ## Load metabolism
    metabolism_fpath = PurePath(asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_fpath)
    
    ## Load fold rules
    fold_rules_db = pd.read_pickle(PurePath("data", "rn2fold", "_db.pkl"))
    fold_rules_path = fold_rules_db.iloc[4]["OUTPUT_PATH"]
    rn2rules = pd.read_pickle(fold_rules_path)
    fold_independent_rns = set()
    foldnet = nf.GlobalFoldNetwork(rn2rules, fold_independent_rns)

    ## Inititalize fold metabolism
    fm = nf.FoldMetabolism(metabolism, foldnet)
    fm.seed_cpds = set((pd.read_csv("data/josh/seed_set.csv")["ID"]))
    fm.seed_folds = set(['spontaneous'])

    current, iteration_dict = fm.rule_order()