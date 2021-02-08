import networkExpansionPy.lib as ne
from networkExpansionPy.folds import FoldRules,fold_expansion
import numpy as np
from random import sample
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--simulation_id", type=str,
                    help="unique id for each random simulation; will be saved to a table on the scc")

# parse the arguments

#outFilePath = args.output

root_path = '/projectnb/bioinfor/SEGRE/goldford/network_expansion/networkExpansionPy'
outFilePath = root_path + '/fold_expansion_scc/results/'
sid = args.simulation_id
foldOrderTableFile = outFilePath + 'fold_table/ft_' + sid + '.csv'
cpdOrderTable = outFilePath + 'cpd_table/ct_' + sid + '.csv'
rxnOrderTable = outFilePath + 'rxn_table/rt_' + sid + '.csv'

# construct global metabolism
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


# define fold rules
fold_rules = FoldRules()
fold_rules.setRules()


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
        

# remove folds that are not included in metabolic network at all
fold_rules.removeFolds(fold_remove)

# define a random permutation on the folds
folds = list(fold_rules.folds)
idx = np.random.permutation(len(folds))
fold_order = [folds[i] for i in idx]


fold_df = pd.DataFrame({'iteration': list(range(len(fold_order))), 'fold': fold_order})
fold_df.to_csv(foldOrderTableFile)


fold_set = set();
cpd_set = seed_set
rxn_set = set();
iteration = 0;
cpds_iteration = {'cid': list(seed_set), 'iteration' : [iteration for x in seed_set]}
rxns_iteration = {'rn': list(rxn_set) , 'iteration' : [iteration for x in rxn_set]}


for fold in fold_order:
    
    iteration = iteration + 1;
    fold_set = fold_set.union(set([fold]))
    c,re,rf = fold_expansion(metabolism,fold_rules,fold_set,cpd_set,rxn_set)
    c = set(c); re = set(re); rf = set(rf)
    
    c_new = [x for x in c if x not in cpd_set]
    re_new = [x for x in re if x not in rxn_set]
    
    c_iter = [iteration for x in c_new]
    r_iter = [iteration for x in re_new]
    cpds_iteration['cid'] = cpds_iteration['cid'] + c_new
    cpds_iteration['iteration'] = cpds_iteration['iteration'] + c_iter
    rxns_iteration['rn'] = rxns_iteration['rn'] + re_new
    rxns_iteration['iteration'] = rxns_iteration['iteration'] + r_iter
    
    # now add seed sets and reaction sets for next iteration
    cpd_set = set(cpd_set).union(c_new)
    rxn_set = rxn_set.union(re_new)

cpd_df = pd.DataFrame(cpds_iteration)
rxn_df = pd.DataFrame(rxns_iteration)
#fold_df = pd.DataFrame({'iteration': list(range(len(fold_order))), 'fold': fold_order})

# join cpd table to fold table
cpd_df = cpd_df.set_index('iteration').join(fold_df.set_index('iteration')).reset_index()
rxn_df = rxn_df.set_index('iteration').join(fold_df.set_index('iteration')).reset_index()
# save dataframes to scc
cpd_df.to_csv(cpdOrderTable)
rxn_df.to_csv(rxnOrderTable)