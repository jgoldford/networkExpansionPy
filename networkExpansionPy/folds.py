from operator import index
import numpy as np
import pandas as pd
from pathlib import PurePath, Path
from copy import copy, deepcopy

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

class GlobalFoldNetwork:

    def __init__(self, rn2rules, fold_independent_rns):

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rules2rn = self.create_foldrules2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rules2rn.keys() for i in fs]) ## all folds
        print("GlobalFoldNetwork initialized\n%i folds available in RUN",len(folds))
        self.fold_independent_rns = fold_independent_rns
        # self.cpds = None ## compounds in fold network only. this would require the mapping from reactions to compounds via metabolism


    def create_foldrules2rn(self, rn2fold):
        fold2rn = dict()
        for rn, rules in rn2fold.items():
            for fs in rules:
                if fs in fold2rn:
                    fold2rn[fs].add(rn)
                else:
                    fold2rn[fs] = set([rn])
        self.rules2rn = fold2rn


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

        self.current_folds = None ## Maybe want to get rid of this
        self.current_cpds = None

        # calculate_scope(cpds_set)

        # self.all_rns = None ## All possible reactions ever; doens't depend on seed set
        ## A few different things
        ## metabolism.rns -- all reactions in kegg universe
        ## fold.rns_all -- all reactions associated with fold universe
        ## fold.rns_i -- all reactions associated with current fold expansion state
        # self.all_cpds = None ## All possible compounds ever; doesnt depend on seed set

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
            self._seed_cpds = seed_cpds 
            self.scope_cpds, self.scope_rns = calculate_scope(self._seed_cpds)
            self.scope_rn2rules = {k:v for k,v in self._f.rn2rules if k in self.scope_rns}
            self.scope_rules2rn = create_foldrules2rn(self.scope_rn2rules)
            self.scope_folds = set([i for fs in self.scope_rules2rn.keys() for i in fs])
            print("Folds in RUN reduced to overlap with scope from fold reactions and fold independent reactions\n%i folds available in RUN",len(folds)) ## scope only includes reactions in folds, and reacitons explicitly input as independent of folds

        else:
            pass 
                
    def calculate_scope(self, seed_cpds):
        ## This is the scope with all reactions in KEGG. 
        ##   But if I wanted I could restrict this to reactions which are only in the folds, or foldindependentreactions, which is a subset of KEGG
        scope_cpds, scope_rns = (set(i) for i in self._m.expand(seed_cpds, reaction_mask=(self._f.rns | self._f.fold_independent_rns)))

    # def filter_to_rules_in_maxiumum_scope(self):
    #     rules2rn = dict()
    #     for k,v in self._f.rules2rn.items():
    #         v_intersect = v & self.scope_rns
    #         if len(v_intersect) > 0:
    #             rules2rn[k] = v_intersect
        
    #     return rules2rn
    def folds2rules(self, folds, rules):
        return {k:v for k,v in rules.items() if k <= rules}
            
    def create_foldrules2rn(self, rn2fold):
        fold2rn = dict()
        for rn, rules in rn2fold.items():
            for fs in rules:
                if fs in fold2rn:
                    fold2rn[fs].add(rn)
                else:
                    fold2rn[fs] = set([rn])
        return fold2rn

    def create_foldrule2subset(fold2rn):
        """
        Given a foldrule -> reaction mapping, identify all foldrules whose corresponding reactions
            are subsets of other foldrule -> reaction mappings.

        Returns mappings which are:
        - subsets (including equal sets)
        - strict subsets (subset must be smaller than equal set)
        - equal sets
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

    def filter_next_iter_to_folds_enabling_new_reactions(self, current_folds):
        """
        Create a mapping of fold:{rns}, where {rns} are the reactions enabled by adding the fold to the existing foldset.

        Mapping discards folds which don't enable any new reactions.

        Used to narrow down possible folds for the next iteration of fold expansion.
        """
        ## I don't think we need to provide existing reactions, because that should be inferabble from the exstiing folds, right?
        ## Well, not all those reactions are reachable from the current state of compounds...

        ## Don't ever try a fold which can't give you any reactions

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
        print("Folds enabling no new reactions during the next NEXT ITERATION removed\n%i folds available for the NEXT ITERATION",len(folds_enabling_no_new_reactions))

        ## need to ouput if a fold doesn't do anything either--that is, by its addition to the existing foldset, it won't be able to enable any new reactions
        return {k:v for k,v in future_fold2rns if k not in folds_enabling_no_new_reactions}

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

    # def find_maximum_scope_folds(self):
    #     for k,v in self.rules2rn.items():
    #         pass

        ## 10/29/2022 I think we want to find a restricted rules2rn mapping that gets rid of reactions unobtainable 
        ##            from the current scope--because unless you're randomly injecting compounds down the road, those folds
        ##            won't actually ever be able to realize those extra reactions
        ## Then, when we do all the narrowing-down in filter_next_iter_to_folds_enabling_new_reactions, we don't want to compare to all_folds
        ##            and rules2rn, we want to compare to the scope-restricted versions of those things.
    
    def next_iter_possible_folds(self, current_folds):

        ## Need to run these two calls every iteration of the fold expansion
        future_fold2rns = self.filter_next_iter_to_folds_enabling_new_reactions(current_folds)
        filtered_folds_to_expand = self.filter_next_iter_to_foldsupersets(future_fold2rns)
        print("Folds whose rules correspond to reactions which are subsets of one another in the next NEXT ITERATION removed\n%i folds available for the NEXT ITERATION",len(filtered_folds_to_expand))
        return filtered_folds_to_expand

    # def folds2rns(self):
    #     """
    #     Get list of reactions possible with current folds (fold independent reactions not included)
    #     """
    #     return set([v for k,v in self.scope_rules2rn.items() if k in self.current_folds])

    # def fold_expansion(self):
    #     current_rns = (self.folds2rn() | self._f.fold_independent_rns)
    #     rxn_set = set(self._m.rxns2tuple(current_rns))
    #     cx,rx = metabolism.expand(self.current_cpds,reaction_mask=rxn_set)
    #     return cx, rx, current_rns

    def fold_expand(self, metabolism, folds, rules2rn, cpds):
        """Doesn't use self"""
        rn_tup_set = set(metabolism.rxns2tuple(set(rules2rns.values())))
        cx,rx = metabolism.expand(cpds,reaction_mask=rn_tup_set)
        return cx, rx


    def effect_per_fold(self, fold, current_folds, current_cpds):
        """
        I think this needs to be modified since i don't have self.current_compounds right now

        self.folds2rules
        self.scope_rules2rn
        self._m
        """

        potential_fold_set = (current_folds | set([fold]))
        potential_rules2rn = self.folds2rules(potential_fold_set, self.scope_rules2rn)
        cx,rx = fold_expand(self._m, potential_fold_set, potential_rules2rn, current_cpds)

        return potential_rules2rn, set(cx), set(rx)

        # old_rules_enabled = self.folds2rules(self.current_folds, self.scope_rules2rn)
        # new_rules_enabled = set(potential_rules2rn.keys()) - set(old_rules_enabled.keys())
        # new_rns_enabled = set(rx) - set([i for s in old_rules_enabled.values() for i in s])
        # new_cpds_enabled = set(cx) - self.current_cpds

        # return new_rules_enabled, new_rns_enabled, new_cpds_enabled

    def loop_through_folds(self, remaining_folds, current_folds, current_cpds):
        """
        Doesn't use self
        """
        next_iter_possible_folds = next_iter_possible_folds(remaining_folds)
        f_effects = dict()
        for f in next_iter_possible_folds:
            _fdict = dict()
            _fdict["rules2rn"], _fdict["cpds"], _fdict["rns"] = effect_per_fold(f, current_folds, current_cpds)
            f_effects[f] = _fdict
        return f_effects

    # define a function that retuns the fold that leads to the largest number of reactions
    # def maxFold(fdict):
    #     q = {fold: len(j) for fold,j in fdict.items()}
    #     q = pd.Series(q)
    #     return q[q == q.max()].index.tolist()

    def maxreactions(self, f_effects):
        k_vcount = {k:len(v["rns"]) for k,v in f_effects.items()}
        return max(k_vcount, key = k_vcount.get)

    def select_next_fold(self, remaining_folds, current_folds, current_cpds, fselect_func=maxreactions):
        f_effects = loop_through_folds(remaining_folds, current_folds, current_cpds)
        next_fold = fselect_func(f_effects)
        return next_fold, f_effects[next_fold]

    ## Need to make sure I run the expansion once without folds first
    ## Need to make sure that the initial "current_folds" includes the seed_folds. 
    ## Can i just store current_folds in my same current dict?

    # def update_current(self, next_fold, fdata, current):

    #     ## Update folds
    #     self.current_folds = (self.current_folds | set([next_fold]))
        
    #     ## Update cpds, rns, rules2rn
    #     current["cpds"] = (fdata["cpds"] | current["cpds"]) 
    #     current["rns"] = (fdata["rns"] | current["rns"]) ## I only need to track these to see what iteration they appeared
    #     current["rules2rn"] = self.folds2rules(self.current_folds, self.scope_rules2rn)

    #     return current


    def fold_order(self):
        
        ## Initial current values
        current = dict() ## current cpds, rns, rules2rn, folds
        current["folds"] = copy.deepcopy(self.seed_folds)
        current["cpds"] = copy.deepcopy(self.seed_cpds)
        current["rules2rn"] = self.folds2rules(current["folds"], self.scope_rules2rn)
        current["rns"] = (self._f.fold_independent_rns | set([i for s in rules2rn_0.values() for i in s]))  ## As written, these are the reactions allowed by the folds, NOT the reactions utilized by the available compounds
        
        ## Iteration data
        keepgoing = True 
        iteration = 0
        iteration_dict = dict()
        iteration_dict["cpds"] = {i:iteration for i in current["cpds"]}
        iteration_dict["rns"] = {i:iteration for i in current["rns"]}
        remaining_folds = (self.scope_folds - current["folds"])
        while keepgoing:
            next_fold, fdata = select_next_fold(remaining_folds, current["folds"], current["cpds"])
            
            ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
            current["folds"] = (current["folds"] | set([next_fold]))
            current["cpds"] = copy.deepcopy(fdata["cpds"])
            current["rules2rn"] = copy.deepcopy(fdata["rules2rn"])
            current["rns"] = copy.deepcopy(fdata["rns"])
            
            ## Store when cpds and rns appear in the expansion
            iteration_dict["cpds"] = {i:iteration for i in current["cpds"] if i not in iteration_dict["cpds"]}
            iteration_dict["rns"] = {i:iteration for i in current["rns"] if i not in iteration_dict["rns"]}

def example_main():
    ## Load metabolism
    metabolism_fpath = PurePath(asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_fpath)
    
    ## Load fold rules
    fold_rules_db = pd.read_pickle(PurePath("data", "rn2fold", "_db.pkl"))
    fold_rules_path = fold_rules_db.iloc[4]["OUTPUT_PATH"]
    rn2rules = pd.read_pickle(fold_rules_path)
    foldnet = GlobalFoldNetwork(rn2rules)

    ## Inititalize fold metabolism
    fm = FoldMetabolism(metabolism, rn2rules)
    fm.seed_cpds = set((pd.read_csv("data/josh/seed_set.csv")["ID"]))
    fm.fold_independent_rns = set()
    fm.seed_folds = set(['spontaneous'])

    fm.expand(foldSet, seedSet, foldIndpendentReactions)
















    ## Define fold rules
    # fold_rules = nf.FoldRules()
    # ecode2rn_path = PurePath("ecode","ecod2rn.keggHMM.17July2021.csv")
    # fold_rules.setRules(path=ecode2rn_path)
    # folds = [list(x) for x in fold_rules.rules.fold_sets.tolist()]
    # folds = [item for sublist in folds for item in sublist]

    # seed_set = list(pd.read_csv("data/josh/seed_set.csv")["ID"])

    # print('starting optimal path finding...')

    # # initialize the algorithm, create a list of folds that remain to be selected during expansion
    # # define folds and reactions that are always available (fold Naive reactions)

    # # define the set of reactions and folds that are accessible initially
    # foldIndpendentReactions = []
    # fold_set = set(['spontaneous']);
    # fold_set_i = set(['spontaneous'])
    # fold_order = {'iteration': [], 'fold': []}
    # # define all the remaining folds
    # folds_remain = set([f for f in folds if f not in fold_set])
    # iteration = 0;

    # # find out all the initially feasible reactions
    # foldEnabledReactions = fold_rules.folds2reactions(fold_set)
    # rxn_set = foldEnabledReactions + foldIndpendentReactions

    # # convert to reaction list to list of tuples with directions
    # rxn_set = set(metabolism.rxns2tuple(rxn_set))

    # # define the initial reaction set that is feasible
    # cpds_iteration = {'cid': list(seed_set), 'iteration' : [iteration for x in seed_set]}
    # rxns_iteration = {'rn': list(rxn_set) , 'iteration' : [iteration for x in rxn_set]}
    # condition = True
    # cpd_set = seed_set

    # # perform network expansion using the restricted set of folds
    # c,re,rf = fold_expansion(metabolism,fold_rules,fold_set_i,cpd_set,rxn_set)


    # foldmetabolism = FoldMetabolism()

#function that peforms fold expansion
def fold_expansion(metabolism,foldRules,fold_set,cpds_set,rxns_seed):
    rn_list = foldRules.folds2reactions(fold_set)
    rn_list = metabolism.rxns2tuple(rn_list) + list(rxns_seed)
    cx,rx = metabolism.expand(cpds_set,reaction_mask=rn_list)
    return cx,rx,rn_list