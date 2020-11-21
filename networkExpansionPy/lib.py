from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
import ray
from random import sample
import os

# define asset path
asset_path,filename = os.path.split(os.path.abspath(__file__))
asset_path = asset_path + '/assets'

def netExp(R,P,x,b):
    k = np.sum(x);
    k0 = 0;
    y = [];
    
    while k > k0:
        k0 = np.sum(x);
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
    return x,y

def netExp_trace(R,P,x,b):
    
    X = []
    Y = []
    
    X.append(x)
    k = np.sum(x);
    k0 = 0;
    y = [];
    
    Y.append(y)
    
    while k > k0:
        k0 = np.sum(x);
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
        X.append(x)
        Y.append(y) 
    return X,Y


def parse_reaction_trace(reaction_trace,network):
    rxns_list = []
    for i in range(1,len(reaction_trace)):
        idx = reaction_trace[i].nonzero()[0]
        rxns = list(network.iloc[:,idx])
        rxns = pd.DataFrame(rxns,columns = ['rn','direction'])
        rxns['iter'] = i
        rxns_list.append(rxns)    
    rxns_list = pd.concat(rxns_list,axis=0)
    return rxns_list

class GlobalMetabolicNetwork:
    
    def __init__(self):
        # load the data
        network = pd.read_csv(asset_path + '/KEGG/network_full.csv')
        cpds = pd.read_csv(asset_path +'/compounds/cpds.tab',sep='\t')
        thermo = pd.read_csv(asset_path +'/reaction_free_energy/kegg_reactions_CC_ph7.0.csv',sep=',')
        self.network = network
        self.compounds = cpds
        self.thermo = thermo
        self.temperature = 50
        self.seedSet = None;
        
    def set_ph(self,pH):
        if ~(type(pH) == str):
            pH = str(pH)
        try:
            thermo = pd.read_csv(asset_path + '/reaction_free_energy/kegg_reactions_CC_ph' + pH + '.csv',sep=',')
            self.thermo = thermo
        except Exception as error:
            print('Failed to open pH files (please use 5.0-9.0 in 0.5 increments)')    
    
    
    def pruneInconsistentReactions(self):
        # remove reactions with qualitatively different sets of elements in reactions and products
        consistent = pd.read_csv(asset_path + '/reaction_sets/reactions_consistent.csv')
        self.network = self.network[self.network.rn.isin(consistent.rn.tolist())]
        
    def pruneUnbalancedReactions(self):
        # only keep reactions that are elementally balanced
        balanced = pd.read_csv(asset_path + '/reaction_sets/reactions_balanced.csv')
        self.network = self.network[self.network.rn.isin(balanced.rn.tolist())]
        
    
    def convertToIrreversible(self):
        network = self.network
        rn = network[['rn']]
        cid = network[['cid']]
        s = network[['s']]
        rn_f = rn + ''
        rn_b = rn + ''
        rn_f['direction'] = 'forward'
        rn_b['direction'] = 'reverse'
        
        nf = cid.join(rn_f).join(s)
        nb = cid.join(rn_b).join(-s)
        self.network = pd.concat([nf,nb],axis=0) 
    
    def setMetaboliteBounds(self,ub = 1e-1,lb = 1e-6): 
        
        self.network['ub'] = ub
        self.network['lb'] = lb
      
    
    def pruneThermodynamicallyInfeasibleReactions(self,keepnan = False):
        fixed_mets = ['C00001','C00080']

        RT = 0.008309424 * (273.15+self.temperature)
        rns  = []
        dirs = []
        dgs = []
        for (rn,direction), dff in self.network.groupby(['rn','direction']):
            effective_deltaG = np.nan
            if rn in self.thermo['!MiriamID::urn:miriam:kegg.reaction'].tolist():
                deltaG = self.thermo[self.thermo['!MiriamID::urn:miriam:kegg.reaction'] == rn]['!dG0_prime (kJ/mol)'].values[0]
                if direction == 'reverse':
                    deltaG = -1*deltaG

                dff = dff[~dff['cid'].isin(fixed_mets)]
                subs = dff[dff['s'] < 0]
                prods = dff[dff['s'] > 0];
                k = np.dot(subs['ub'].apply(np.log),subs['s']) + np.dot(prods['lb'].apply(np.log),prods['s'])

                effective_deltaG = RT*k + deltaG

            dgs.append(effective_deltaG)
            dirs.append(direction)
            rns.append(rn)

        res = pd.DataFrame({'rn':rns,'direction':dirs,'effDeltaG':dgs})
        if ~keepnan:
            res = res.dropna()
        
        res = res[res['effDeltaG'] < 0].set_index(['rn','direction'])
        res = res.drop('effDeltaG',axis=1)
        self.network = res.join(self.network.set_index(['rn','direction'])).reset_index()
    
    def initialize_metabolite_vector(self,seedSet):
        if self.seedSet is None:
            print('No seed set')
        else:
            network = self.network.pivot_table(index='cid',columns = ['rn','direction'],values='s').fillna(0)
            x0 = np.array([x in seedSet for x in network.index.get_level_values(0)]) * 1;        
            return x0
        
    

