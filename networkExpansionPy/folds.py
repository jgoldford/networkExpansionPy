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

    def __init__(self, rn2rules):

        self.rn2rules = rn2rules ## all rns to fold rules
        self.rules2rn = self.create_foldrules2rn(rn2rules) ## all fold rules to rns
        self.rns = set(rn2rules.keys()) ## reactions in fold network only
        self.folds = set([i for fs in self.rules2rn.keys() for i in fs]) ## all folds
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

    def __init__(self, metabolism, rn2rules, preexpansion=False):
        # load the data
        self.metabolism = metabolism ## Metabolism object
        self.rn2rules = rn2rules ## all rns to fold rules
        self.rules2rn = self.create_foldrules2rn(rn2rules) ## all fold rules to rns
        self.all_folds = set([i for fs in self.rules2rn.keys() for i in fs]) ## All possible folds ever; doesn't depent on seed set
        self.current_folds = None

        self.seed_folds = None
        self.seed_cpds = None
        self.maximum_scope_calculated = False
        self.maximum_scope_cpds = None
        self.maximum_scope_rns = None
        self.maximum_scope_folds = None

        # maximum_scope(cpds_set)

        # self.all_rns = None ## All possible reactions ever; doens't depend on seed set
        ## A few different things
        ## metabolism.rns -- all reactions in kegg universe
        ## fold.rns_all -- all reactions associated with fold universe
        ## fold.rns_i -- all reactions associated with current fold expansion state
        # self.all_cpds = None ## All possible compounds ever; doesnt depend on seed set
        

    def maximum_scope(self, seed_cpds):
        self.maximum_scope_cpds, self.maximum_scope_rns = self.metabolism.expand(seed_cpds)
        self.maximum_scope_calculated = True

    def create_foldrules2rn(self, rn2fold):
        fold2rn = dict()
        for rn, rules in rn2fold.items():
            for fs in rules:
                if fs in fold2rn:
                    fold2rn[fs].add(rn)
                else:
                    fold2rn[fs] = set([rn])
        self.rules2rn = fold2rn

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

    def find_superset_folds(self, future_fold2rns):
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

    def find_folds_enabling_new_reactions(self):
        """
        Create a mapping of fold:{rns}, where {rns} are the reactions enabled by adding the fold to the existing foldset

        Mapping discards folds which don't enable any new reactions
        """
        ## I don't think we need to provide existing reactions, because that should be inferabble from the exstiing folds, right?
        ## Well, not all those reactions are reachable from the current state of compounds...

        ## Don't ever try a fold which can't give you any reactions

        ## Don't ever try a fold which can only give you a subset of reactions that adding a different fold could give you, unless you care about ties, in which case you should flag that as a kwarg.
        # fold2foldsubset, fold2foldsmallersubset = create_foldrule2subset(rules2rn)
        
        future_fold2rns = dict()
        for fold in self.all_folds-self.current_folds:
            future_folds = self.current_folds | set([fold])
            future_foldrules2rns = {k:v for k,v in self.rules2rn.items() if k <= future_folds}
            future_rns = set([rn for v in future_foldrules2rns.values() for rn in v])
            future_fold2rns[fold] = future_rns

        current_fold2rn = {k:v for k,v in self.rules2rn.items() if k <= self.current_folds}
        current_rns = set([rn for v in current_fold2rn.values() for rn in v])
        folds_enabling_no_new_reactions = set([k for k,v in future_fold2rns.items() if v==current_rns])

        ## need to ouput if a fold doesn't do anything either--that is, by its addition to the existing foldset, it won't be able to enable any new reactions
        return {k:v for k,v in future_fold2rns if k not in folds_enabling_no_new_reactions}

    def find_maximum_scope_folds(self):
        for k,v in self.rules2rn.items():
            pass

        ## 10/29/2022 I think we want to find a restricted rules2rn mapping that gets rid of reactions unobtainable 
        ##            from the current scope--because unless you're randomly injecting compounds down the road, those folds
        ##            won't actually ever be able to realize those extra reactions
        ## Then, when we do all the narrowing-down in find_folds_enabling_new_reactions, we don't want to compare to all_folds
        ##            and rules2rn, we want to compare to the scope-restricted versions of those things.
    
    def prepare_for_expansion(self):
        folds_to_expand = self.find_superset_folds(self.find_folds_enabling_new_reactions())


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
    fm.seed_folds = 
















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