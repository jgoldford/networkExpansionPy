import networkExpansionPy.lib as ne
from networkExpansionPy.folds import FoldRules,fold_expansion
import numpy as np
from random import sample
import pandas as pd
import argparse

import warnings
from scipy.sparse import (spdiags, SparseEfficiencyWarning, csc_matrix,
    csr_matrix, isspmatrix, dok_matrix, lil_matrix, bsr_matrix)
warnings.simplefilter('ignore',SparseEfficiencyWarning)

# define some helper functions
def maxFold(fdict):
    q = {fold: len(j) for fold,j in fdict.items()}
    q = pd.Series(q)
    return q[q == q.max()].index.tolist()


root_path = '/projectnb/bioinfor/SEGRE/goldford/network_expansion/networkExpansionPy'
#root_path = '/Users/joshuagoldford/Documents/github/networkExpansionPy'

outFilePath = root_path + '/fold_expansion_scc/results/'

foldOrderTable = outFilePath + 'optimal_path_ec4/foldOrderTable.csv'
cpdOrderTable = outFilePath + 'optimal_path_ec4/cpdOrderTable.csv'
rxnOrderTable = outFilePath + 'optimal_path_ec4/rnOrderTable.csv'

foldOrderTableTmp = outFilePath + 'optimal_path_ec4/foldOrderTable.tmp.csv'
cpdOrderTableTmp = outFilePath + 'optimal_path_ec4/cpdOrderTable.tmp.csv'
rxnOrderTableTmp = outFilePath + 'optimal_path_ec4/rnOrderTable.tmp.csv'
# construct global metabolism

print('building metabolic network...')
metabolism = ne.GlobalMetabolicNetwork()
metabolism.convertToIrreversible()

# remove all O2 dependent reactions
oxygen_dependent_rxns = metabolism.network[metabolism.network.cid.isin(['C00007'])].rn.unique().tolist()
o2_independent_rxns = [x for x in metabolism.network.rn.unique().tolist() if x not in oxygen_dependent_rxns]

# only keep anaerobic reactions
metabolism.subnetwork(o2_independent_rxns)



# define seed compounds
cpds = pd.read_csv(root_path + '/networkExpansionPy/assets/compounds/seeds.csv')
cpds['CID'] = cpds['CID'].apply(lambda x: x.strip())
seed_set = cpds['CID'].unique().tolist()


print('building fold-to-reaction mapping...')
# define fold rules
fold_rules = FoldRules()
fold_rules.setRules(path='ecode/ecod2rn.ec4.12Feb2021.csv')


# determine scope of folds (reaction sets that involve the fold)
fold_dict = {}
for fold in fold_rules.folds:
    foldset = set([fold])
    rns = fold_rules.rules[fold_rules.rules.fold_sets.apply(lambda x: foldset.issubset(x))].rn.unique().tolist()
    fold_dict[fold] = rns
    
    
# if the fold cover reactions that are not in the metabolic netowrk, do not use fold in iterative fold expansion code
mrxns = metabolism.network.rn.unique().tolist()
fold_remove = []
for key,values in fold_dict.items():
    if all([x not in mrxns for x in values]):
        fold_remove.append(key)

# remove all rules with folds that are not used in metabolic network or are erroneous
fold_remove = fold_remove + ['PDBChainNotFound']
# remove folds that are not included in metabolic network at all
fold_rules.removeFolds(fold_remove)


# define fold list
folds = list(fold_rules.folds)
folds = [f for f in folds if f.lower() != 'pdbchainnotfound']
#folds = ['2002','304','2007','2003','2004','230']


fold_set = set();
cpd_set = seed_set
rxn_set = set();
iteration = 0;
cpds_iteration = {'cid': list(seed_set), 'iteration' : [iteration for x in seed_set]}
rxns_iteration = {'rn': list(rxn_set) , 'iteration' : [iteration for x in rxn_set]}

print('starting optimal path finding...')

# initialize the algorithm with a fold set remain
fold_order = {'iteration': [], 'fold': []}
folds_remain = folds
condition = True


while condition:
    iteration = iteration + 1;
    cpd_dict = {}
    rn_dict = {}

    # for each reaction fold in the fold remain step, run network expansion with the fold
    for fold in folds_remain:
        
        fold_set_i = fold_set.union(set([fold]))
        # run fold expansion
        c,re,rf = fold_expansion(metabolism,fold_rules,fold_set_i,cpd_set,rxn_set)
        c = set(c); re = set(re); rf = set(rf)
        
        # save new reactions and metabolites to dictionary
        c_new = [x for x in c if x not in cpd_set]
        re_new = [x for x in re if x not in rxn_set]
        cpd_dict[fold] = c_new
        rn_dict[fold] = re_new
        print('finished with fold: ' + fold + ', R: ' + str(len(re_new)) + ' C: ' +  str(len(c_new)))


    
    max_fold_set = maxFold(rn_dict)

    if len(max_fold_set) > 1:
        cpd_dict_subset = {k: cpd_dict[k] for k in max_fold_set}
        max_fold_set = maxFold(cpd_dict_subset)
        if len(max_fold_set) > 1:
            max_fold = sample(max_fold_set,1)[0]
        else:
            max_fold = max_fold_set[0]
    else:
        max_fold = max_fold_set[0]

    cpd_new = set(cpd_dict[max_fold])
    rn_new = set(rn_dict[max_fold])

    c_iter = [iteration for x in cpd_new]
    r_iter = [iteration for x in rn_new]
    cpds_iteration['cid'] = cpds_iteration['cid'] + list(cpd_new)
    cpds_iteration['iteration'] = cpds_iteration['iteration'] + c_iter
    rxns_iteration['rn'] = rxns_iteration['rn'] + list(rn_new)
    rxns_iteration['iteration'] = rxns_iteration['iteration'] + r_iter

    # now add seed sets and reaction sets for next iteration
    cpd_set = set(cpd_set).union(cpd_new)
    rxn_set = rxn_set.union(rn_new)
    fold_set = fold_set.union(set([max_fold]))

    fold_order['fold'].append(max_fold)
    fold_order['iteration'].append(iteration)

    # update the fold remain list
    folds_remain = [fold for fold in folds_remain if fold not in fold_order['fold']]
    # stop the algorithm if (a) largest network is of size 0, or (b) no more folds are possible
    if (len(cpd_new) + len(rn_new)) < 1:
        condition = False
    if len(folds_remain) < 1:
        condition = False
    print('finished with iteration: ' + str(iteration))

    cpd_df_tmp = pd.DataFrame(cpds_iteration)
    rxn_df_tmp = pd.DataFrame(rxns_iteration)
    fold_df_tmp = pd.DataFrame(fold_order)
    cpd_df_tmp.to_csv(cpdOrderTableTmp)
    rxn_df_tmp.to_csv(rxnOrderTableTmp)
    fold_df_tmp.to_csv(foldOrderTableTmp)


cpd_df = pd.DataFrame(cpds_iteration)
rxn_df = pd.DataFrame(rxns_iteration)
fold_df = pd.DataFrame(fold_order)

# save dataframes to scc
cpd_df.to_csv(cpdOrderTable)
rxn_df.to_csv(rxnOrderTable)
fold_df.to_csv(foldOrderTable)

#cpd_df.to_csv('res.cpd.csv')
#rxn_df.to_csv('res.rxn.csv')
#fold_df.to_csv('res.foldOrder.csv')