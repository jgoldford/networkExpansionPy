from operator import index
import numpy as np
import pandas as pd
from pathlib import PurePath, Path
from copy import copy, deepcopy
import timeit
from pprint import pprint

asset_path = PurePath(__file__).parent / "assets"

# define a function that determins if a reaction rule (x) is feasible with foldSet
def rule2rn(foldSet,x):
    if x.issubset(foldSet):
    #if foldSet.issubset(x)
        return True
    else:
        return False

# define a function that returns a list of reactions that are feasible with foldSet.  rules_sub contains all the fold set rules for all reactions
def folds2rn(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: rule2rn(foldSet,x))
    y = rules_sub[feasible]
    if len(y) > 0:
        rns = rules_sub[feasible].rn.unique().tolist()
    else:
        rns = []
    return rns

# same as above but return only the set of rules that are feasible
def folds2rules(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: rule2rn(foldSet,x))
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
                
        rns = folds2rn(self.rules,foldSet)
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

class GlobalFoldNetwork:

    def __init__(self, rn2rules, fold_independent_rns):

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rules2rn = self.create_foldrules2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rules2rn.keys() for i in fs]) ## all folds
        print("GlobalFoldNetwork initialized\n%i folds available in RUN"%(len(self.folds)))
        self.fold_independent_rns = fold_independent_rns

    def create_foldrules2rn(self, rn2fold):
        fold2rn = dict()
        for rn, rules in rn2fold.items():
            for fs in rules:
                if fs in fold2rn:
                    fold2rn[fs].add(rn)
                else:
                    fold2rn[fs] = set([rn])
        return fold2rn


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
            self.scope_rules2rn = self.create_foldrules2rn(self.scope_rn2rules)
            self.scope_folds = set([i for fs in self.scope_rules2rn.keys() for i in fs])
            print("...done.")
            print("Folds in RUN reduced to overlap with scope from fold reactions and fold independent reactions\n%i folds available in RUN"%len(self.scope_folds)) ## scope only includes reactions in folds, and reacitons explicitly input as independent of folds

        else:
            pass 
                
    def calculate_scope(self, seed_cpds):
        """
        Calculate the scope of the seed_cpds with all reactions enabled by the global fold network (including fold_independent reactions)
        """
        rn_tup_set = set(self._m.rxns2tuple((self._f.rns | self._f.fold_independent_rns)))
        scope_cpds, scope_rns = self._m.expand(seed_cpds, reaction_mask=rn_tup_set)
        return set(scope_cpds), set([i[0] for i in scope_rns])

    def folds2rules(self, folds, rules):
        """Doesn't use self"""
        return {k:v for k,v in rules.items() if k <= set(folds)}
            
    def create_foldrules2rn(self, rn2fold):
        """Doesn't use self"""
        fold2rn = dict()
        for rn, rules in rn2fold.items():
            for fs in rules:
                if fs in fold2rn:
                    fold2rn[fs].add(rn)
                else:
                    fold2rn[fs] = set([rn])
        return fold2rn

    def create_foldrule2subset(self,fold2rn):
        """
        Given a foldrule -> reaction mapping, identify all foldrules whose corresponding reactions
            are subsets of other foldrule -> reaction mappings.

        Returns mappings which are:
        - subsets (including equal sets)
        - strict subsets (subset must be smaller than equal set)
        - equal sets

        Doesn't use self
        """
        fold2foldsubset = dict()
        fold2foldstrictsubset = dict()
        fold2equalfold = dict()
        for k, v in fold2rn.items():
            fold2foldsubset[k] = set()
            fold2foldstrictsubset[k] = set()
            fold2equalfold[k] = set()
            for k2, v2 in fold2rn.items():
                if k != k2:
                    if v2 <= v:
                        fold2foldsubset[k].add(k2)
                        if v2!=v:
                            fold2foldstrictsubset[k].add(k2)
                        else:
                            fold2equalfold[k].add(k2)
        
        return fold2foldsubset, fold2foldstrictsubset, fold2equalfold
    ##############################
    def filter_next_iter_to_rules_enabling_new_reactions(self, current_folds):
        current_fold2rn = {k:v for k,v in self.scope_rules2rn.items() if k <= current_folds}
        current_rns = set([rn for v in current_fold2rn.values() for rn in v])
        # return {k:v for k,v in self.scope_rules2rn.items() if not v <= current_rns} ## rules2rn only which enable new reactions
        ## to make this in the same format as below, need k=rule, v=reactions enabled by that rule plus all reactions which are already enabled.
        ## e.g.
        return {k:(v | current_rns) for k,v in self.scope_rules2rn.items() if not v <= current_rns}

    def create_equal_rule_groups(self, fold2equalfold):

        # create set of equal fold rule groups, will be a set of frozensets of frozensets (!). Middle frozenset contains equivilent rules. Inner most frozenset contains a single rule.
        equal_rule_groups = []
        for k,v in fold2equalfold.items():
            equal_rule_group = [k]
            for fs in v:
                if len(fs)>0:
                    equal_rule_group.append(fs)
            equal_rule_groups.append(frozenset(equal_rule_group))        
        equal_rule_groups = set(equal_rule_groups)
        ## go back through and make everything into sorted tuples for repeatability
        return [sorted([sorted(i) for i in group]) for group in equal_rule_groups]

    # def all_strictsuperset_rules(self, fold2foldstrictsubset):
    #     strictsubset_rules = set([frozenset(j) for i in fold2foldstrictsubset.values() for j in i])
    #     return set([i for i in fold2foldstrictsubset.keys() if i not in strictsubset_rules])
    
    def filter_next_iter_to_rulesupersets(self, future_fold2rns):
        _, fold2foldstrictsubset, fold2equalfold = self.create_foldrule2subset(future_fold2rns)
        
        equal_rule_groups = self.create_equal_rule_groups(fold2equalfold)

        strictsubset_rules = set([frozenset(j) for i in fold2foldstrictsubset.values() for j in i])
        strictsuperset_rules = set([i for i in fold2foldstrictsubset.keys() if i not in strictsubset_rules]) 
        
        # print("strictsubset_rules: ", strictsubset_rules)
        # print("equal_rule_groups: ", equal_rule_groups)
        ## This is is a sorted list of lists, no sets involved anymore. But strictsuperset_rules is still a set of frozensets
        equal_rule_groups_that_arent_subsets = [i for i in equal_rule_groups if set(i[0]) not in strictsubset_rules] # only need to check 1 fold per group since they're equal
        
        ## I should modify the below line to always return the smallest len() equivilent rule. BUT it is possible that a longer rule will include just 1 required new fold, so it should be chosen above a 2-fold rule where both folds are new.
        ## Maybe need to break this into two functions, one returns the strict supersets, one returns equivilent rules
        ## for equivilent rules in groups, choose rule where (rule_fold_i - current_folds) are shortest in length. if tie, sort rules, and sort groups, and choose first
        # smallest_number_fold_from_each_group = {i[0] for i in equal_rule_groups_that_arent_subsets} # replicatable sorting

        ## Only need to expand 1. strict superset folds, and 2. one fold from each equal_fold_group, provided that group is not in strictsubset_folds
        return strictsuperset_rules, equal_rule_groups_that_arent_subsets

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

    def remove_current_folds_from_rules(self, current_folds, rules):
        # if type(rules)!=set:
        # print("current folds: ",current_folds)
        # print("rules: ",rules)
        new_rules = [frozenset(i) - current_folds for i in rules]
        return [i for i in new_rules if len(i)>0]

    def remove_current_folds_from_equal_rule_groups(self, current_folds, rule_groups):
        # if type(rules)!=set:
        # print("current folds: ",current_folds)
        # print("rule_groups: ",rule_groups)
        new_rule_groups = []
        for rg in rule_groups:
            new_group = []
            for i in rg:
                new_rule = frozenset(i) - current_folds
                if len(i) > 0:
                    new_group.append(new_rule)
            if len(new_group) > 0:
                new_rule_groups.append(new_group)
        print(new_rule_groups)
        return new_rule_groups
    
    def next_iter_possible_rules(self, current_folds):

        ## Need to run these two calls every iteration of the fold expansion
        future_fold2rns = self.filter_next_iter_to_rules_enabling_new_reactions(current_folds)
        strictsuperset_rules, equal_rule_groups_that_arent_subsets = self.filter_next_iter_to_rulesupersets(future_fold2rns)
        ## Filter to strictsuperset_rules that only require a single additional fold from the current set
        ##      Only try other ones if single fold additions don't yeild any new reactions
        ## Filter to equal_rule_groups_that_arent_subsets that only require a single additional fold from the current set
        ##      Only try other ones if single fold additions don't yeild any new reactions
        strictsuperset_rule_dict = self.organize_set_elements_by_size(self.remove_current_folds_from_rules(current_folds, strictsuperset_rules))
        equal_rule_dict = self.organize_equal_rule_groups_by_size(self.remove_current_folds_from_equal_rule_groups(current_folds, equal_rule_groups_that_arent_subsets))
        print(equal_rule_dict)
        ## How do i handle making the equal_rule dict? i think after I create the current_fold_from_equal_rule_groups
        ##      I need to make a dict within each equal_rule_group. Done.
        ## equal_rule_dict is now a list of dicts, where each dict contains equivilent rule sets

        # print("-> Folds whose rules correspond to reactions which are subsets of one another in the next NEXT ITERATION removed\n-> ... %i folds available for the NEXT ITERATION"%len(filtered_folds_to_expand))
        return strictsuperset_rule_dict, equal_rule_dict #filtered_folds_to_expand

    def effect_per_rule(self, rule, current_folds, current_cpds):
        """
        self.folds2rules
        self.scope_rules2rn
        self._m
        """

        potential_fold_set = (current_folds | set([rule]))
        # print(f"{potential_fold_set=}")
        potential_rules2rn = self.folds2rules(potential_fold_set, self.scope_rules2rn)
        # print(f"{potential_rules2rn=}")
        cx,rx = self.fold_expand(self._m, potential_fold_set, potential_rules2rn, self._f.fold_independent_rns, current_cpds)

        return potential_rules2rn, cx, rx #set(cx), set(rx)

    def loop_through_rules(self, current_folds, current_cpds, current_rns):

        strictsuperset_rule_dict, equal_rule_dict = self.next_iter_possible_rules(current_folds)

        rule_sizes = set(list(strictsuperset_rule_dict.keys()) + [i for d in equal_rule_dict for i in d.keys()])

        for fsize in sorted(rule_sizes):
            f_effects = dict()
            if fsize in strictsuperset_rule_dict:
                for rule in strictsuperset_rule_dict[fsize]:
                    _fdict = dict()
                    _fdict["rules2rn"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule(rule, current_folds, current_cpds)
                    # print(f"{f_effects=}")
                    if len(_fdict["rns"] - current_rns) > 0:
                        f_effects[f] = _fdict
            for d in equal_rule_dict:
                if fsize in d:
                    for rule in d[fsize]:
                        _fdict = dict()
                        _fdict["rules2rn"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule(rule, current_folds, current_cpds)
                        # print(f"{f_effects=}")
                        if len(_fdict["rns"] - current_rns) > 0:
                            f_effects[f] = _fdict

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(f_effects) > 1:
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
        # print("feffects:")
        # pprint(f_effects)
        next_rule = fselect_func(f_effects)
        return next_rule, f_effects[next_rule]

    def rule_order(self):

        if (self.seed_cpds == None) or (self.seed_folds == None):
            raise ValueError("self.seed_cpds and self.seed_folds must not be None")
        
        ## Initialize current values
        current = {
            "folds": deepcopy(self.seed_folds),
            "cpds": deepcopy(self.seed_cpds),
            "rns":set([])
            }
        
        ## Initialize iteration data
        keepgoing = True 
        iteration = 0
        iteration_dict = {
            "cpds":dict(), 
            "rns":dict(), 
            "folds":{"fold_independent":0}
            }
        
        ## Avoid updating folds on the 0th iteration since they don't apply until iteration=1
        iteration_dict = self.update_iteration_dict(iteration_dict, {k:v for k,v in current.items() if k!="folds"}, iteration)
        remaining_folds = (self.scope_folds - current["folds"])
        iteration+=1

        # print("iteration dict INIT: ", iteration_dict)
        # print("current dict INIT: ", current)

        ## First expansion (using only seed folds and fold independent reactions)
        init_rules2rn = self.folds2rules(current["folds"], self.scope_rules2rn)
        current["cpds"], current["rns"] = self.fold_expand(self._m, current["folds"], init_rules2rn, self._f.fold_independent_rns, current["cpds"])
        iteration_dict = self.update_iteration_dict(iteration_dict, current, iteration)

        # print("iteration dict after first expansion: ", iteration_dict)
        # print("current dict after first expansion: ", current)


        while keepgoing:
            start = timeit.default_timer()
            iteration += 1
            print("\nITERATION: ", iteration)
            print("maximum n folds remaining: ", len(remaining_folds))
            print("maximum remaining_folds:\n", remaining_folds)
            next_rule, fdata = self.select_next_rule(current["folds"], current["cpds"], current["rns"])
            print("next fold: ", next_rule)
            print("fdata:")
            pprint(fdata)
            exec_time = timeit.default_timer() - start
            print("iteration runtime: ", exec_time)
            remaining_folds = (remaining_folds - set(next_rule))
            if len(remaining_folds) == 0:
                keepgoing = False

            ## Stop conditions
            if (fdata["cpds"] == current["cpds"]) and (fdata["rns"] == current["rns"]):
                keepgoing = False    

                # ## If no reactions were added, try multiple folds simultaneously
                # multifold_rules2rn = remaining_multifold_rules()     
                # if len(multifold_rules2rn) > 0:
                #     ## Need to basically do the existing process to select fold, except selecting all folds of an entire rule (of length 2) instead of selecting just 1 fold
                #     ##      Always take smallest number of folds which provide any new reactions/compounds. It might also be the case that after adding a double fold, single folds again
                #     ##      become viable.
                #     next_rule, fdata = self.select_next_rule(current["folds"], current["cpds"])

            else:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                current["folds"] = (current["folds"] | set(next_rule))
                current["cpds"] = fdata["cpds"]
                current["rns"] = fdata["rns"]
                
                ## Store when cpds and rns appear in the expansion
                iteration_dict = self.update_iteration_dict(iteration_dict, current, iteration)

        return current, iteration_dict

    ##############################        

    def filter_next_iter_to_folds_enabling_new_reactions(self, current_folds):
        """
        Create a mapping of fold:{rns}, where {rns} are the reactions enabled by adding the fold to the existing foldset.

        Mapping discards folds which don't enable any new reactions.

        Used to narrow down possible folds for the next iteration of fold expansion.
        """
        ## I don't think we need to provide existing reactions, because that should be inferabble from the exstiing folds, right?
        ## Well, not all those reactions are reachable from the current state of compounds...

        ## Don't ever try a fold which can only give you a subset of reactions that adding a different fold could give you, unless you care about ties, in which case you should flag that as a kwarg.
        # fold2foldsubset, fold2foldsmallersubset = create_foldrule2subset(rules2rn)
        
        future_fold2rns = dict()
        for fold in self.scope_folds-current_folds:
            future_folds = current_folds | set([fold])
            future_foldrules2rns = {k:v for k,v in self.scope_rules2rn.items() if k <= future_folds}
            future_rns = set([rn for v in future_foldrules2rns.values() for rn in v])
            future_fold2rns[fold] = future_rns

        current_fold2rn = {k:v for k,v in self.scope_rules2rn.items() if k <= current_folds}
        current_rns = set([rn for v in current_fold2rn.values() for rn in v])
        folds_enabling_no_new_reactions = set([k for k,v in future_fold2rns.items() if v==current_rns])
        print("-> Folds enabling no new reactions during the next NEXT ITERATION removed\n-> ... %i folds available for the NEXT ITERATION"%len(set(future_fold2rns.keys()) - folds_enabling_no_new_reactions))

        ## need to ouput if a fold doesn't do anything either--that is, by its addition to the existing foldset, it won't be able to enable any new reactions
        return {k:v for k,v in future_fold2rns.items() if k not in folds_enabling_no_new_reactions}

    def filter_next_iter_to_foldsupersets(self, future_fold2rns):
        """
        Narrow down which future folds need to be expanded based on redundancy of reactions enabled across folds.

        Returns set of folds.

        :param future_fold2rns: dictionary of fold:{rns}, where {rns} are the reactions enabled by adding the fold to the existing foldset
        """

        _, fold2foldstrictsubset, fold2equalfold = self.create_foldrule2subset(future_fold2rns)

        equal_fold_groups = {tuple({k}.union(v)) for k,v in fold2equalfold.items() if len(v)>0} # create set of equal fold groups
        equal_fold_groups = {tuple(sorted(int(i) for i in g)) for g in equal_fold_groups} # turn folds into ints for sorting
        equal_fold_groups = {tuple(str(i) for i in g) for g in equal_fold_groups} # turn back into strs

        strictsubset_folds = set([j for i in fold2foldstrictsubset.values() for j in i])
        strictsuperset_folds = set([i for i in fold2foldstrictsubset.keys() if i not in strictsubset_folds])
        
        equal_fold_groups_that_arent_subsets = {i for i in equal_fold_groups if i[0] not in strictsubset_folds} # only need to check 1 fold per group since they're equal
        smallest_number_fold_from_each_group = {i[0] for i in equal_fold_groups_that_arent_subsets} # replicatable sorting

        ## Only need to expand 1. strict superset folds, and 2. one fold from each equal_fold_group, provided that group is not in strictsubset_folds
        return strictsuperset_folds | smallest_number_fold_from_each_group
    
    def next_iter_possible_folds(self, current_folds):

        ## Need to run these two calls every iteration of the fold expansion
        future_fold2rns = self.filter_next_iter_to_folds_enabling_new_reactions(current_folds)
        filtered_folds_to_expand = self.filter_next_iter_to_foldsupersets(future_fold2rns)
        print("-> Folds whose rules correspond to reactions which are subsets of one another in the next NEXT ITERATION removed\n-> ... %i folds available for the NEXT ITERATION"%len(filtered_folds_to_expand))
        return filtered_folds_to_expand

    def fold_expand(self, metabolism, folds, rules2rn, fold_independent_rns, cpds):
        """Doesn't use self"""
        fold_rns = set([i for s in rules2rn.values() for i in s])
        rn_tup_set = set(metabolism.rxns2tuple(fold_rns | fold_independent_rns))
        cx,rx = metabolism.expand(cpds,reaction_mask=rn_tup_set)
        return set(cx), set([i[0] for i in rx])


    def effect_per_fold(self, fold, current_folds, current_cpds):
        """
        self.folds2rules
        self.scope_rules2rn
        self._m
        """

        potential_fold_set = (current_folds | set([fold]))
        # print(f"{potential_fold_set=}")
        potential_rules2rn = self.folds2rules(potential_fold_set, self.scope_rules2rn)
        # print(f"{potential_rules2rn=}")
        cx,rx = self.fold_expand(self._m, potential_fold_set, potential_rules2rn, self._f.fold_independent_rns, current_cpds)

        return potential_rules2rn, cx, rx #set(cx), set(rx)

    def loop_through_folds(self, current_folds, current_cpds):
        next_iter_possible_folds = self.next_iter_possible_folds(current_folds)
        f_effects = dict()
        for f in next_iter_possible_folds:
            _fdict = dict()
            _fdict["rules2rn"], _fdict["cpds"], _fdict["rns"] = self.effect_per_fold(f, current_folds, current_cpds)
            f_effects[f] = _fdict
        print(f"{f_effects=}")
        return f_effects

    # def maxreactions(f_effects):
    #     """Doesn't use self"""
    #     k_vcount = {k:len(v["rns"]) for k,v in f_effects.items()}
    #     k_vcount = dict(sorted(k_vcount.items())) ## Sort for reproduceability
    #     return max(k_vcount, key = k_vcount.get)

    def select_next_fold(self, current_folds, current_cpds, fselect_func=maxreactions):
        """Doesn't use self"""
        f_effects = self.loop_through_folds(current_folds, current_cpds)
        # print("feffects:")
        # pprint(f_effects)
        next_fold = fselect_func(f_effects)
        return next_fold, f_effects[next_fold]

    def update_iteration_dict(self, iteration_dict, current, iteration):
        for dtype, ids in current.items():
            for i in ids:
                if i not in iteration_dict[dtype]:
                    iteration_dict[dtype][i] = iteration
        return iteration_dict

    def fold_order(self):

        if (self.seed_cpds == None) or (self.seed_folds == None):
            raise ValueError("self.seed_cpds and self.seed_folds must not be None")
        
        ## Initialize current values
        current = {
            "folds": deepcopy(self.seed_folds),
            "cpds": deepcopy(self.seed_cpds),
            "rns":set([])
            }
        
        ## Initialize iteration data
        keepgoing = True 
        iteration = 0
        iteration_dict = {
            "cpds":dict(), 
            "rns":dict(), 
            "folds":{"fold_independent":0}
            }
        
        ## Avoid updating folds on the 0th iteration since they don't apply until iteration=1
        iteration_dict = self.update_iteration_dict(iteration_dict, {k:v for k,v in current.items() if k!="folds"}, iteration)
        remaining_folds = (self.scope_folds - current["folds"])
        iteration+=1

        print("iteration dict INIT: ", iteration_dict)
        print("current dict INIT: ", current)

        ## First expansion (using only seed folds and fold independent reactions)
        init_rules2rn = self.folds2rules(current["folds"], self.scope_rules2rn)
        current["cpds"], current["rns"] = self.fold_expand(self._m, current["folds"], init_rules2rn, self._f.fold_independent_rns, current["cpds"])
        iteration_dict = self.update_iteration_dict(iteration_dict, current, iteration)

        print("iteration dict after first expansion: ", iteration_dict)
        print("current dict after first expansion: ", current)


        while keepgoing:
            start = timeit.default_timer()
            iteration += 1
            print("\nITERATION: ", iteration)
            print("maximum n folds remaining: ", len(remaining_folds))
            print("maximum remaining_folds:\n", remaining_folds)
            next_fold, fdata = self.select_next_fold(current["folds"], current["cpds"])
            print("next fold: ", next_fold)
            print("fdata:")
            pprint(fdata)
            exec_time = timeit.default_timer() - start
            print("iteration runtime: ", exec_time)
            remaining_folds = (remaining_folds - set([next_fold]))
            if len(remaining_folds) == 0:
                keepgoing = False

            ## Stop conditions
            if (fdata["cpds"] == current["cpds"]) and (fdata["rns"] == current["rns"]):
                keepgoing = False    

                ## If no reactions were added, try multiple folds simultaneously
                multifold_rules2rn = remaining_multifold_rules()     
                if len(multifold_rules2rn) > 0:
                    ## Need to basically do the existing process to select fold, except selecting all folds of an entire rule (of length 2) instead of selecting just 1 fold
                    ##      Always take smallest number of folds which provide any new reactions/compounds. It might also be the case that after adding a double fold, single folds again
                    ##      become viable.
                    next_fold, fdata = self.select_next_fold(current["folds"], current["cpds"])

            else:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                current["folds"] = (current["folds"] | set([next_fold]))
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

    current, iteration_dict = fm.fold_order()