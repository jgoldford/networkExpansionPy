import numpy as np
import pandas as pd
import os
from copy import copy, deepcopy

# define asset path
asset_path,filename = os.path.split(os.path.abspath(__file__))
asset_path = asset_path + '/assets'


def rule2rn(foldSet,x):
    if x.issubset(foldSet):
    #if foldSet.issubset(x)
        return True
    else:
        return False
    
def folds2rn(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: rule2rn(foldSet,x))
    y = rules_sub[feasible]
    if len(y) > 0:
        rns = rules_sub[feasible].rn.unique().tolist()
    else:
        rns = []
    return rns
    
def folds2rules(rules_sub,foldSet):
    feasible = rules_sub['fold_sets'].apply(lambda x: rule2rn(foldSet,x))
    rule_df = rules_sub[feasible]
    return rule_df


class FoldRules:
    
    def __init__(self):
        # load the data
        self.rules = None;
        self.rn = None
        self.folds = None;
        
    def copy(self):
        return deepcopy(self)
        
    def setRules(self,path = '/ecode/ecod2rn.ec3.07Feb2021.csv'):
        rules = pd.read_csv(asset_path + path)
        self.rns = rules.rn.unique().tolist()
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
def fold_expansion(metabolism,foldRules,fold_set,cpd_set,rxns_seed):
    rxns_feasible = foldRules.folds2reactions(fold_set)
    rxns_total = list(rxns_feasible) + list(rxns_seed)
    m = metabolism.copy()
    m.subnetwork(rxns_total)
    cpds_ne,rxns_ne = m.expand(list(cpd_set))
    return cpds_ne,rxns_ne,rxns_feasible