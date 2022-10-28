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
        self.allrules = setRules(path=path) ## This should be all rules, every step of expansion will subset
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

# define a function to run network expansion using a fold set
# depends on metabolism object and foldRules objects
# OLD FOLD EXPANSION CODE: REPLACED 10/5 W FASTER VERSION
#def fold_expansion(metabolism,foldRules,fold_set,cpd_set,rxns_seed):
#    rxns_feasible = foldRules.folds2reactions(fold_set)
#    rxns_total = list(rxns_feasible) + list(rxns_seed)
#    m = metabolism.copy()
#    m.subnetwork(rxns_total)
#    cpds_ne,rxns_ne = m.expand(list(cpd_set))
#    return cpds_ne,rxns_ne,rxns_feasible

class FoldMetabolism:

    def __init__(self, metabolism, foldRules, preexpansion=False):
        # load the data
        self.metabolism = metabolism
        self.foldRules = foldRules ## This could just be my pickled dict of rns->folds
        self.seed_folds = seed_folds
        self.maximum_scope_calculated = False
        self.maximum_scope_cpds = set()
        self.maximum_scope_rns = set()

        maximum_scope(cpds_set)

    def maximum_scope(self, cpds_set):
        self.maximum_scope_cpds, self.maximum_scope_rns = self.metabolism.expand(cpds_set)
        self.maximum_scope_calculated = True


def example_main():
    ## Define metabolism
    metabolism_fpath = PurePath(ne.asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_fpath)
    
    ## Define fold rules
    fold_rules = nf.FoldRules()
    ecode2rn_path = PurePath("ecode","ecod2rn.keggHMM.17July2021.csv")
    fold_rules.setRules(path=ecode2rn_path)
    folds = [list(x) for x in fold_rules.rules.fold_sets.tolist()]
    folds = [item for sublist in folds for item in sublist]

    seed_set = list(pd.read_csv("data/josh/seed_set.csv")["ID"])

    print('starting optimal path finding...')

    # initialize the algorithm, create a list of folds that remain to be selected during expansion
    # define folds and reactions that are always available (fold Naive reactions)

    # define the set of reactions and folds that are accessible initially
    foldIndpendentReactions = []
    fold_set = set(['spontaneous']);
    fold_set_i = set(['spontaneous'])
    fold_order = {'iteration': [], 'fold': []}
    # define all the remaining folds
    folds_remain = set([f for f in folds if f not in fold_set])
    iteration = 0;

    # find out all the initially feasible reactions
    foldEnabledReactions = fold_rules.folds2reactions(fold_set)
    rxn_set = foldEnabledReactions + foldIndpendentReactions

    # convert to reaction list to list of tuples with directions
    rxn_set = set(metabolism.rxns2tuple(rxn_set))

    # define the initial reaction set that is feasible
    cpds_iteration = {'cid': list(seed_set), 'iteration' : [iteration for x in seed_set]}
    rxns_iteration = {'rn': list(rxn_set) , 'iteration' : [iteration for x in rxn_set]}
    condition = True
    cpd_set = seed_set

    # perform network expansion using the restricted set of folds
    c,re,rf = fold_expansion(metabolism,fold_rules,fold_set_i,cpd_set,rxn_set)


    foldmetabolism = FoldMetabolism()

#function that peforms fold expansion
def fold_expansion(metabolism,foldRules,fold_set,cpds_set,rxns_seed):
    rn_list = foldRules.folds2reactions(fold_set)
    rn_list = metabolism.rxns2tuple(rn_list) + list(rxns_seed)
    cx,rx = metabolism.expand(cpds_set,reaction_mask=rn_list)
    return cx,rx,rn_list