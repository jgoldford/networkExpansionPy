from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
from random import sample
import os
import json
from copy import copy, deepcopy
import zipfile
import pickle

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



def netContract(R,P,b,x_ac,y,yext):
    # code for running network contraction algorithm
    Xac = []
    Yac = []
    
    #x_ac = x.copy()
    Xac.append(x_ac)
    k0 = np.sum(x_ac);
    n_reactions = np.size(R,1)
    #yactive = csr_matrix(np.multiply((y.toarray()),1-yext.toarray()))   
    yactive = yext.__lt__(1).multiply(y)
    k0 = yactive.sum()
    Yac.append(yactive)
    k = 0
    while k < k0:
        # find
        k0 = yactive.sum();
        x_ex = np.dot(P,yext).astype('bool').astype('int')
        x_ac = np.dot(P,yactive).astype('bool').astype('int')
        # find all extinct metabolites
        #x_ex = csr_matrix(np.multiply((x_ex.toarray()),1-x_ac.toarray()))
        
        # metabolite has to be not active and extinct
        x_ex = x_ac.__lt__(1).multiply(x_ex)
        
        #compute extinct reactions
        yext = np.dot(R.transpose(),x_ex).astype('bool').astype('int')
        # compute active reactions
        
        #yactive = csr_matrix(np.multiply((yactive.toarray()),1-yext.toarray()))   
        yactive = yext.__lt__(1).multiply(yactive)
        #k = np.sum(x_ac);
        k = yactive.sum()
        Xac.append(x_ac)
        Yac.append(yactive)
        
        
    return Xac,Yac


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

def load_json_network(rdict):
    network_list = []
    consistent_rids = []
    for rid,v in rdict.items():
        if v["glycans"] == False: ## This skips all reactions with glycans
            cids = v["left"] + v["right"]
            try: ## This skips all reactions with n stoichiometries
                stoichs = [-int(i) for i in v["left_stoichiometries"]]+[int(i) for i in v["right_stoichiometries"]]
                network_list+=list(zip(cids,[rid for _ in range(len(stoichs))],stoichs))
            except:
                pass
            if v["element_conservation"]==True:
                consistent_rids.append(rid)
    return pd.DataFrame(network_list,columns=("cid","rn","s")), pd.DataFrame(consistent_rids,columns=["rn"])

def _load_tuple_network(tlist):
    """
    Load a simple network, defined by 2-tuples of reactant lists and product lists
    
    Note: Intended for DEV only
    Note: Stoichiometry information in output is only accurate for directionality
    
        tlist_example = [
             (["A","B"],["C"]),
             (["C","D"],["E","F"]),
             (["E","F"],["G"]),
             (["G","H"],["I"]),
             (["A","J"],["I"])]
    """
    rows = list()
    for i,d in enumerate(tlist):
        if len(d)!=2: raise ValueError("reactions must be two-tuples")
        for cid in d[0]:
            rows.append({"rn":i, "cid":cid, "s":-1})
        for cid in d[1]:
            rows.append({"rn":i, "cid":cid, "s":1})
    return pd.DataFrame(rows)

class GlobalMetabolicNetwork:
    
    def __init__(self,metabolism="KEGG_OG"):
        # load the data        
        if metabolism == "KEGG_OG":
            network = pd.read_csv(asset_path + '/KEGG/network_full.csv')
            cpds = pd.read_csv(asset_path +'/compounds/cpds.txt',sep='\t')
            thermo = pd.read_csv(asset_path +'/reaction_free_energy/kegg_reactions_CC_ph7.0.csv',sep=',')
            self.network = network
            self.thermo = thermo
            self.compounds = cpds ## Includes many compounds without reactions
            
        
        elif metabolism == "ecg":
            with open(os.path.join(asset_path,"ecg","master_from_kegg_2021-01-05.json")) as f:
                ecg = json.load(f)
            network, consistent_rxns = load_ecg_network(ecg)
            self.network = network
            self.consistent_rxns = consistent_rxns
            self.compounds = pd.DataFrame(self.network["cid"].unique(),columns=["cid"]) ## Only includes compounds with reactions

        elif metabolism == "KEGG":
            with zipfile.ZipFile(os.path.join(asset_path,"KEGG","2021.05.31-18.06.52","entries_detailed","reaction.json.zip"),"r") as z:
                rdict = json.loads(z.read(z.infolist()[0]).decode())
            network, consistent_rxns = load_json_network(rdict)
            self.network = network
            self.consistent_rxns = consistent_rxns
            self.compounds = pd.DataFrame(self.network["cid"].unique(),columns=["cid"]) ## Only includes compounds with reactions

        elif metabolism == "dev":
            ## Just for testing, etc.
            self.network = None

        else:
            raise(ValueError("'metabolism' must be one of 'KEGG_OG, 'ecg', 'KEGG'"))

        self.metabolism = metabolism
        self.temperature = 25
        self.seedSet = None
        self.rid_to_idx = None
        self.idx_to_rid = None
        self.cid_to_idx = None
        self.idx_to_cid = None
        self.S = None
        
    def copy(self):
        return deepcopy(self)
        
    def set_ph(self,pH):
        if ~(type(pH) == str):
            pH = str(pH)
        if self.metabolism == "KEGG_OG":
            try:
                thermo = pd.read_csv(asset_path + '/reaction_free_energy/kegg_reactions_CC_ph' + pH + '.csv',sep=',')
                self.thermo = thermo
            except Exception as error:
                print('Failed to open pH files (please use 5.0-9.0 in 0.5 increments)')    
        elif self.metabolism == "ecg":
            try:
                self.thermo = self.load_ecg_thermo(self.metabolism,pH)
            except:
                raise ValueError("Try another pH, that one appears not to be in the ecg json")
        else:
            raise(NotImplementedError("pH not yet implemented for metabolism = %s"%self.metabolism)) 


    def load_ecg_thermo(self,ph=9):
        thermo_list = []
        for rid,v in self.metabolism["reactions"].items():
            
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

    
    def pruneInconsistentReactions(self):
        # remove reactions with qualitatively different sets of elements in reactions and products
        if self.metabolism=="KEGG_OG":
            consistent = pd.read_csv(asset_path + '/reaction_sets/reactions_consistent.csv')
            self.network = self.network[self.network.rn.isin(consistent.rn.tolist())]
        elif self.metabolism=="KEGG" or self.metabolism=="ecg":
            self.network = self.network[self.network.rn.isin(self.consistent_rxns.rn.tolist())]
        else:
            raise(NotImplementedError("Function not yet implemented for metabolism = %s"%self.metabolism)) 

    def pruneUnbalancedReactions(self):
        # only keep reactions that are elementally balanced
        if self.metabolism=="KEGG_OG":
            balanced = pd.read_csv(asset_path + '/reaction_sets/reactions_balanced.csv')
            self.network = self.network[self.network.rn.isin(balanced.rn.tolist())]
        else:
            raise(NotImplementedError("Function not yet implemented for metabolism = %s"%self.metabolism)) 
        
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

        if not hasattr(self, 'thermo'):
            raise(AttributeError("Metabolism has no thermo data."))

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
        if keepnan:
            # change effective free energy to negative number, so it passes the next filter
            res['effDeltaG'] = res['effDeltaG'].fillna(-1)
        else:
            res = res.dropna()

        #res = res[res['effDeltaG'] < 0].set_index(['rn','direction'])
        res = res[~(res['effDeltaG'] > 0)].set_index(['rn','direction'])
        res = res.drop('effDeltaG',axis=1)
        self.network = res.join(self.network.set_index(['rn','direction'])).reset_index()
    
    def pruneReactionsFromMetabolite(self,cpds):
        # find all reactions that use metabolites in cpds list
        reactions_to_drop = self.network[self.network.cid.isin(cpds)].rn.unique().tolist()
        reactions_to_keep = [x for x in self.network.rn.unique().tolist() if x not in reactions_to_drop]
        # only keep reactions that do not use that metabolite
        self.subnetwork(reactions_to_keep)

    def initialize_metabolite_vector(self,seedSet):
        if seedSet is None:
            print('No seed set')
        else:
            x0 = np.zeros([len(self.cid_to_idx)],dtype=int)
            for x in set(seedSet)&set(self.cid_to_idx.keys()):
                x0[self.cid_to_idx[x]] = 1     
            return x0

    def initialize_reaction_vector(self,reactionSet):
        if reactionSet is None:
            print('No reactions in set')
        else:
            x0 = np.zeros([len(self.rid_to_idx)],dtype=int)
            for x in set(reactionSet)&set(self.rid_to_idx.keys()):
                x0[self.rid_to_idx[x]] = 1     
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


    def create_RP_from_irreversible_network(self):
        
        # only use entries in teh network dataframe with s < 0:
        reactant_network = self.network[self.network.s<0]

        R = np.zeros([len(self.cid_to_idx),len(self.rid_to_idx)])
        
        for c,r,d,s in zip(reactant_network["cid"],reactant_network["rn"],reactant_network["direction"],reactant_network["s"]):
            R[self.cid_to_idx[c],self.rid_to_idx[(r,d)]] = 1

        # only use entries in teh network dataframe with s > 0:
        product_network = self.network[self.network.s>0]

        P = np.zeros([len(self.cid_to_idx),len(self.rid_to_idx)])
        
        for c,r,d,s in zip(product_network["cid"],product_network["rn"],product_network["direction"],product_network["s"]):
            P[self.cid_to_idx[c],self.rid_to_idx[(r,d)]] = 1

        return R,P

    def create_iteration_dict(self,M,idx_to_id):
        idx_iter = dict()
        for i,row in enumerate(M):
            idxs = np.nonzero(row.toarray().T[0])[0]
            for idx in idxs:
                if idx not in idx_iter:
                    idx_iter[idx] = i

        id_iter = dict()
        for idx,i in idx_iter.items():
            id_iter[idx_to_id[idx]] = i  

        return id_iter
        
    def expand(self,seedSet,algorithm='naive',reaction_mask = None):
        # constructre network from skinny table and create matricies for NE algorithm
        # if (self.rid_to_idx is None) or (self.idx_to_rid is None):
        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        # if (self.cid_to_idx is None) or (self.idx_to_cid is None):
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        # if self.S is None:
        #self.S = self.create_S_from_irreversible_network()
        x0 = self.initialize_metabolite_vector(seedSet)
        
        #R = (self.S < 0)*1
        #P = (self.S > 0)*1
        R,P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        # add a new term that uses sparse matrix multiplication for R and P to zero out reactions are that are not accessible
        if reaction_mask is not None:
            reaction_mask = self.initialize_reaction_vector(reaction_mask)
            reaction_mask = csr_matrix(np.diag(reaction_mask))
            P = P*reaction_mask
            R = R*reaction_mask

        x0 = csr_matrix(x0)
        x0 = x0.transpose()
        if algorithm.lower() == 'naive':
            x,y = netExp(R,P,x0,b)
        elif algorithm.lower() == 'cr':
            x,y = netExp_cr(R,P,x0,b)
        elif algorithm.lower() == 'trace':
            X,Y = netExp_trace(R,P,x0,b)
        else:
            raise ValueError('algorithm needs to be naive (compound stopping criteria) or cr (reaction/compound stopping criteria)')
        
        if algorithm.lower() == 'trace':
    
            compound_iteration_dict = self.create_iteration_dict(X,self.idx_to_cid)
            reaction_iteration_dict = self.create_iteration_dict(Y,self.idx_to_rid)
            return compound_iteration_dict, reaction_iteration_dict

        else:
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


    def contract(self,seedSet,reactionScope,compoundScope,extinctReactions):
        # constructre network from skinny table and create matricies for NE algorithm
        # if (self.rid_to_idx is None) or (self.idx_to_rid is None):
        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        # if (self.cid_to_idx is None) or (self.idx_to_cid is None):
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        # if self.S is None:
        #self.S = self.create_S_from_irreversible_network()
        
        # create vectors for reactionScope, compoundScope
        xactive = self.initialize_metabolite_vector(compoundScope)
        yactive = self.initialize_reaction_vector(reactionScope)
        yextinct = self.initialize_reaction_vector(extinctReactions)
        #R = (self.S < 0)*1
        #P = (self.S > 0)*1
        R,P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        #x0 = csr_matrix(x0)
        #x0 = x0.transpose()
        xactive = csr_matrix(xactive).transpose()
        yactive = csr_matrix(yactive).transpose()
        yextinct = csr_matrix(yextinct).transpose()
        # reun contraction algorithm
        X,Y = netContract(R,P,b,xactive,yactive,yextinct)
        x = X[-1]
        y = Y[-1]

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

    def run_expansions(self,seedSets):
        # constructre network from skinny table and create matricies for NE algorithm
        # if (self.rid_to_idx is None) or (self.idx_to_rid is None):
        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        # if (self.cid_to_idx is None) or (self.idx_to_cid is None):
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        # if self.S is None:
        #self.S = self.create_S_from_irreversible_network()
        
        #R = (self.S < 0)*1
        #P = (self.S > 0)*1
        R,P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        compoundScopes = []
        reactionScopes = []
        for seedSet in seedSets:
            x0 = self.initialize_metabolite_vector(seedSet)
            x0 = csr_matrix(x0)
            x0 = x0.transpose()
            x,y = netExp(R,P,x0,b)
            
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
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)
            return compoundScopes,reactionScopes

    def run_contractions(self,seedSet,reactionScope,compoundScope,extinctReactionSets):
        # constructre network from skinny table and create matricies for NC algorithm
        # if (self.rid_to_idx is None) or (self.idx_to_rid is None):

        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        # if (self.cid_to_idx is None) or (self.idx_to_cid is None):
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        # if self.S is None:
        #self.S = self.create_S_from_irreversible_network()
        
        # create vectors for reactionScope, compoundScope
        xactive = self.initialize_metabolite_vector(compoundScope)
        yactive = self.initialize_reaction_vector(reactionScope)
        
        #R = (self.S < 0)*1
        #P = (self.S > 0)*1
        R,P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        #x0 = csr_matrix(x0)
        #x0 = x0.transpose()
        xactive = csr_matrix(xactive).transpose()
        yactive = csr_matrix(yactive).transpose()
     

        compoundScopes = []
        reactionScopes = []

        for extinctReactions in extinctReactionSets:
            yextinct = self.initialize_reaction_vector(extinctReactions)
            yextinct = csr_matrix(yextinct).transpose()
            # run contraction algorithm
            X,Y = netContract(R,P,b,xactive,yactive,yextinct)
            x = X[-1]
            y = Y[-1]            
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
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)
        
        return compoundScopes,reactionScopes


    def run_expansions_reactionMasks(self,seedSet,maskedReactionSets):
        # run many expansions, but masking with reaction list in the lists maskedReactionSets
        R,P = self.create_RP_from_irreversible_network()
        b = sum(R)
        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        x0 = self.initialize_metabolite_vector(seedSet)
        x0 = csr_matrix(x0)
        x0 = x0.transpose()

        compoundScopes = []
        reactionScopes = []

        for rxns_removed in maskedReactionSets:
            # build new R and P matricies with Masks
            yextinct = self.initialize_reaction_vector(rxns_removed)
            reaction_mask = csr_matrix(np.diag(1-yextinct))
            Pstar = P*reaction_mask
            Rstar = R*reaction_mask
            # run contraction algorithm
            x,y = netExp(Rstar,Pstar,x0,b)

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
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)
        
        return compoundScopes,reactionScopes


    def ne_output_to_graph(self,cpds,rxns):
        # build a network constructing metabolites from prior iteration to connecting subsequent iteration
        # input: cpds (rxns) is a dict with compound id (reaction id) as the key, and iteration as the value
        # output: dataframe with target and source node and iteration for  

        graph = {'target': [], 'source': [], 'iteration':[]}
        maxIter = np.array(list(cpds.values())).max()
        
        cpd_iter_df = pd.DataFrame({'cid': [x[0] for x in cpds.items()],'iter': [x[1] for x in cpds.items()]})

        for i in range(maxIter,0,-1):
            molecules = [x[0] for x in cpds.items() if x[1] == i]
            reactions = [x[0] for x in rxns.items() if x[1] == i]
            reactions = self.network.set_index(['rn','direction']).loc[reactions]
            for molecule in molecules:
                # reactions
                # find any reactions that produce this metabolite
                r_sample = reactions[ ( reactions.cid == molecule) & (reactions.s>0)].sample(1)
                # find metabolite in reaction that was produced in previous iteration
                r_sample_cpds = reactions.loc[r_sample.index].set_index('cid').join(cpd_iter_df.set_index('cid'))
                r_sample_cpds = r_sample_cpds[r_sample_cpds.iter == i-1].sample(1)
                molecules_origin = r_sample_cpds.index.tolist()[0]
                graph['target'].append(molecule)
                graph['source'].append(molecules_origin)
                graph['iteration'].append(i)

        graph = pd.DataFrame(graph)
        return graph

    def save(self,name):
        path_to_save = asset_path + '/metabolic_networks/' + name + ".pkl"
        with open(path_to_save, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def rxns2tuple(self,rn_list):
        t = self.network[self.network.rn.isin(rn_list)][['rn','direction']].drop_duplicates()
        rn_list_tuple = list(zip(t.rn.tolist(),t.direction.tolist()))
        return rn_list_tuple
    