from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
import ray
from random import sample
import os
import json
from copy import copy, deepcopy

# define asset path
asset_path,filename = os.path.split(os.path.abspath(__file__))
asset_path = asset_path + '/assets'

def netExp(R,P,x,b):
    k = np.sum(x);
    k0 = 0;
    n_reactions = np.size(R,1)
    y = csr_matrix(np.zeros(n_reactions))
    while k > k0:
        k0 = np.sum(x);
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
    return x,y

# define a new network expansion, s.t. stopping criteria is now no new compounds or reactions can be added at subsequent iterations
def netExp_cr(R,P,x,b):
    k = np.sum(x);
    k0 = 0;
    n_reactions = np.size(R,1)
    y = csr_matrix(np.zeros(n_reactions))
    l = 0
    l0 = 0;
        
    while (k > k0) | (l > l0):
        k0 = np.sum(x);
        l0 = np.sum(y)
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
        l = np.sum(y)
    return x,y


def netExp_trace(R,P,x,b):
    
    X = []
    Y = []
    
    X.append(x)
    k = np.sum(x);
    k0 = 0;
    n_reactions = np.size(R,1)
    y = csr_matrix(np.zeros(n_reactions))
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


def isRxnCoenzymeCoupled(rxn,cosubstrate,coproduct):
    g = rxn[rxn.cid.isin([cosubstrate,coproduct])]
    out = False
    if len(g) > 1:
        if g.s.sum() == 0:
            out = True
    return out

def load_ecg_network(ecg):
    network_list = []
    consistent_rids = []
    for rid,v in ecg["reactions"].items():
        cids = v["left"] + v["right"]
        try: ## This skips all reactions with n stoichiometries
            stoichs = [-int(i) for i in v["metadata"]["left_stoichiometries"]]+[int(i) for i in v["metadata"]["right_stoichiometries"]]
            network_list+=list(zip(cids,[rid for _ in range(len(stoichs))],stoichs))
        except:
            pass
        if v["metadata"]["element_conservation"]==True:
            consistent_rids.append(rid)
    return pd.DataFrame(network_list,columns=("cid","rn","s")), pd.DataFrame(consistent_rids,columns=["rn"])

def load_ecg_thermo(ecg,ph=9):
    thermo_list = []
    for rid,v in ecg["reactions"].items():
        
        phkey = str(ph)+"pH_100mM"
        
        if v["metadata"]["dg"][phkey]["standard_dg_prime_value"] == None:
            dg = np.nan
        else:
            dg = v["metadata"]["dg"][phkey]["standard_dg_prime_value"]
            
        if v["metadata"]["dg"][phkey]["standard_dg_prime_error"] == None:
            dgerror = np.nan
        else:
            dgerror = v["metadata"]["dg"][phkey]["standard_dg_prime_error"]
            
        if v["metadata"]["dg"][phkey]["is_uncertain"] == None:
            note = "uncertainty is too high"
        else:
            note = np.nan

        thermo_list.append((rid,
            dg,
            dgerror,
            v["metadata"]["dg"][phkey]["p_h"],
            v["metadata"]["dg"][phkey]["ionic_strength"]/1000,
            v["metadata"]["dg"][phkey]["temperature"],
            note)) 

    return pd.DataFrame(thermo_list, columns = ("!MiriamID::urn:miriam:kegg.reaction","!dG0_prime (kJ/mol)","!sigma[dG0] (kJ/mol)","!pH","!I (mM)","!T (Kelvin)","!Note"))         

class GlobalMetabolicNetwork:
    
    def __init__(self,ecg_json=None):
        # load the data
        if ecg_json == None:
            network = pd.read_csv(asset_path + '/KEGG/network_full.csv')
            cpds = pd.read_csv(asset_path +'/compounds/cpds.txt',sep='\t')
            thermo = pd.read_csv(asset_path +'/reaction_free_energy/kegg_reactions_CC_ph7.0.csv',sep=',')
            self.compounds = cpds
            self.ecg = None
        else:
            with open(ecg_json) as f:
                ecg = json.load(f)
            network, consistent_rxns = load_ecg_network(ecg)
            thermo = load_ecg_thermo(ecg)
            self.ecg = ecg
            self.consistent_rxns = consistent_rxns
            ## self.compounds appears not to be used so it's not included here

        self.network = network
        self.thermo = thermo
        self.temperature = 25
        self.seedSet = None
        self.rid_to_idx = None
        self.idx_to_rid = None
        self.cid_to_idx = None
        self.idx_to_cid = None
        
    def copy(self):
        return deepcopy(self)
        
    def set_ph(self,pH):
        if ~(type(pH) == str):
            pH = str(pH)
        if self.ecg == None:
            try:
                thermo = pd.read_csv(asset_path + '/reaction_free_energy/kegg_reactions_CC_ph' + pH + '.csv',sep=',')
                self.thermo = thermo
            except Exception as error:
                print('Failed to open pH files (please use 5.0-9.0 in 0.5 increments)')    
        else:
            try:
                self.thermo = load_ecg_thermo(self.ecg,pH)
            except:
                raise ValueError("Try another pH, that one appears not to be in the ecg json")

    
    def pruneInconsistentReactions(self):
        # remove reactions with qualitatively different sets of elements in reactions and products
        if self.ecg==None:
            consistent = pd.read_csv(asset_path + '/reaction_sets/reactions_consistent.csv')
            self.network = self.network[self.network.rn.isin(consistent.rn.tolist())]
        else:
            self.network = self.network[self.network.rn.isin(self.consistent_rxns.rn.tolist())]
        
    def pruneUnbalancedReactions(self):
        # only keep reactions that are elementally balanced
        balanced = pd.read_csv(asset_path + '/reaction_sets/reactions_balanced.csv')
        self.network = self.network[self.network.rn.isin(balanced.rn.tolist())]
        
    def subnetwork(self,rxns):
        # only keep reactions that are in list
        self.network = self.network[self.network.rn.isin(rxns)]
        
    def addGenericCoenzymes(self):
        replace_metabolites = {'C00003': 'Generic_oxidant', 'C00004': 'Generic_reductant', 'C00006': 'Generic_oxidant',  'C00005': 'Generic_reductant','C00016': 'Generic_oxidant','C01352':'Generic_reductant'}
        coenzyme_pairs = {}
        coenzyme_pairs['NAD'] = ['C00003','C00004']
        coenzyme_pairs['NADP'] = ['C00006','C00005']
        coenzyme_pairs['FAD'] = ['C00016','C01352']
        coenzyme_pairs = pd.DataFrame(coenzyme_pairs).T.reset_index()
        coenzyme_pairs.columns = ['id','oxidant','reductant']
        # create reactions copies with coenzyme pairs
        new_rxns = []
        new_thermo = [];
        for idx,rxn in self.network.groupby('rn'):
            z = any([isRxnCoenzymeCoupled(rxn,row.oxidant,row.reductant) for x,row in coenzyme_pairs.iterrows()])
            if z:
                new_rxn = rxn.replace(replace_metabolites).groupby(['cid','rn']).sum().reset_index()
                new_rxn['rn'] = new_rxn['rn'] = idx + '_G'
                new_rxns.append(new_rxn)
                t = self.thermo[self.thermo['!MiriamID::urn:miriam:kegg.reaction'] == idx].replace({idx:  idx + '_G'})
                new_thermo.append(t)

        new_rxns = pd.concat(new_rxns,axis=0)
        new_thermo = pd.concat(new_thermo,axis=0)

        self.network = pd.concat([self.network,new_rxns],axis=0)
        self.thermo = pd.concat([self.thermo,new_thermo],axis=0)

    
    def convertToIrreversible(self):
        nf = self.network.copy()
        nb = self.network.copy()
        nf['direction'] = 'forward'
        nb['direction'] = 'reverse'
        nb['s'] = -nb['s']
        net = pd.concat([nf,nb],axis=0)
        net = net.set_index(['cid','rn','direction']).reset_index()
        self.network = net
    
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
        
        #res = res[res['effDeltaG'] < 0].set_index(['rn','direction'])
        res = res[~(res['effDeltaG'] > 0)].set_index(['rn','direction'])
        res = res.drop('effDeltaG',axis=1)
        self.network = res.join(self.network.set_index(['rn','direction'])).reset_index()
    
    def initialize_metabolite_vector(self,seedSet):
        if seedSet is None:
            print('No seed set')
        else:
            x0 = np.zeros([len(self.cid_to_idx)],dtype=int)
            for x in seedSet:
                x0[self.cid_to_idx[x]] = 1     
            return x0

    def create_reaction_dicts(self):
        rids = set(zip(self.network["rn"],self.network["direction"]))
        rid_to_idx = dict()
        idx_to_rid = dict()
        for v, k in enumerate(rids):
            rid_to_idx[k] = v
            idx_to_rid[v] = k
        
        return rid_to_idx, idx_to_rid

    def create_compound_dicts(self):
        cids = set(self.network["cid"])
        cid_to_idx = dict()
        idx_to_cid = dict()
        for v, k in enumerate(cids):
            cid_to_idx[k] = v
            idx_to_cid[v] = k
        
        return cid_to_idx, idx_to_cid

    def create_S_from_irreversible_network(self):
        
        S = np.zeros([len(self.cid_to_idx),len(self.rid_to_idx)])
            
        for c,r,d,s in zip(self.network["cid"],self.network["rn"],self.network["direction"],self.network["s"]):
            S[self.cid_to_idx[c],self.rid_to_idx[(r,d)]] = s

        return S
        
    def expand(self,seedSet,algorithm='naive'):
        # constructre network from skinny table and create matricies for NE algorithm
        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        x0 = self.initialize_metabolite_vector(seedSet)
        S = self.create_S_from_irreversible_network()
        R = (S < 0)*1
        P = (S > 0)*1
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        x0 = csr_matrix(x0)
        x0 = x0.transpose()
        if algorithm.lower() == 'naive':
            x,y = netExp(R,P,x0,b)
        elif algorithm.lower() == 'cr':
            x,y = netExp_cr(R,P,x0,b)
        else:
            raise ValueError('algorithm needs to be naive (compound stopping criteria) or cr (reaction/compound stopping criteria)')
        
        # convert to list of metabolite ids and reaction ids
        if x.toarray().sum() > 0:
            cidx = np.nonzero(x.toarray().T[0])[0]
            compounds = [self.idx_to_cid[i] for i in cidx]
        else:
            compounds = []
            
        if y.toarray().sum() > 0:
            ridx = np.nonzero(y.toarray().T[0])[0]
            reactions = [self.idx_to_rid[i] for i in ridx]
        else:
            reactions = [];
            
        return compounds,reactions

