import unittest
import networkExpansionPy.lib as ne
import pandas as pd
from scipy.sparse import csr_matrix
# from pandas.testing import assert_frame_equal

def dfs_have_equal_content(df1,df2):
    ## Ignores row order and index
    return df1.sort_values(by=list(df1.columns),axis=0).reset_index(drop=True).equals(df2.sort_values(by=list(df2.columns),axis=0).reset_index(drop=True))


class TestGlobalMetabolicNetworkInit(unittest.TestCase):

    def test_kegg_og_no_args(self):
        kegg_og = ne.GlobalMetabolicNetwork()
        self.assertEqual(kegg_og.metabolism,"KEGG_OG")
        self.assertIsInstance(kegg_og.network,pd.core.frame.DataFrame)
        
        nsize = len(kegg_og.network)
        self.assertGreater(nsize,1)
        
        kegg_og.pruneUnbalancedReactions()
        nsize_ubpruned = len(kegg_og.network)
        self.assertLessEqual(nsize_ubpruned,nsize)
        
        kegg_og.pruneInconsistentReactions()
        nsize_inconpruned = len(kegg_og.network)
        self.assertLessEqual(nsize_inconpruned,nsize_ubpruned)

        kegg_og.set_ph(7.0)
        self.assertIsInstance(kegg_og.thermo,pd.core.frame.DataFrame)

        kegg_og.convertToIrreversible()
        self.assertEqual(len(kegg_og.network),2*nsize_inconpruned)

        kegg_og.setMetaboliteBounds(ub=1e-1,lb=1e-6)
        self.assertIn("ub",kegg_og.network)
        self.assertIn("lb",kegg_og.network)

        nsize_prethermoprune = len(kegg_og.network)
        kegg_og.pruneThermodynamicallyInfeasibleReactions(keepnan=False)
        self.assertLessEqual(len(kegg_og.network),nsize_prethermoprune)

    def test_ecg(self):
        ecg = ne.GlobalMetabolicNetwork("ecg")
        self.assertEqual(ecg.metabolism,"ecg")
        self.assertIsInstance(ecg.network,pd.core.frame.DataFrame)
        
        nsize = len(ecg.network)
        self.assertGreater(nsize,1)
        
        with self.assertRaises(NotImplementedError):
            ecg.pruneUnbalancedReactions()

        ecg.pruneInconsistentReactions()
        nsize_inconpruned = len(ecg.network)
        self.assertLessEqual(nsize_inconpruned,nsize)

        with self.assertRaises(ValueError):
            ecg.set_ph(7.0)

        ecg.convertToIrreversible()
        self.assertEqual(len(ecg.network),2*nsize_inconpruned)

        ecg.setMetaboliteBounds(ub=1e-1,lb=1e-6)
        self.assertIn("ub",ecg.network)
        self.assertIn("lb",ecg.network)

        with self.assertRaises(AttributeError):
            ecg.pruneThermodynamicallyInfeasibleReactions()

    def test_KEGG(self):
        kegg = ne.GlobalMetabolicNetwork("KEGG")
        self.assertEqual(kegg.metabolism,"KEGG")
        self.assertIsInstance(kegg.network,pd.core.frame.DataFrame)
        
        nsize = len(kegg.network)
        self.assertGreater(nsize,1)
        
        with self.assertRaises(NotImplementedError):
            kegg.pruneUnbalancedReactions()

        kegg.pruneInconsistentReactions()
        nsize_inconpruned = len(kegg.network)
        self.assertLessEqual(nsize_inconpruned,nsize)

        with self.assertRaises(NotImplementedError):
            kegg.set_ph(7.0)

        kegg.convertToIrreversible()
        self.assertEqual(len(kegg.network),2*nsize_inconpruned)

        kegg.setMetaboliteBounds(ub=1e-1,lb=1e-6)
        self.assertIn("ub",kegg.network)
        self.assertIn("lb",kegg.network)

        with self.assertRaises(AttributeError):
            kegg.pruneThermodynamicallyInfeasibleReactions()

class Test_create_iteration_dict(unittest.TestCase):

    # @classmethod
    # def setup_class(self):
    def setUp(self):
        self.kegg = ne.GlobalMetabolicNetwork("KEGG")
        self.kegg.pruneInconsistentReactions()
        self.kegg.convertToIrreversible()
        self.seedSet = ["C00001",
                "C00011",
                "C00014",
                "C00033",
                "C00058",
                "C00283",
                "C00288",
                "C00697"]

        self.kegg.rid_to_idx, self.kegg.idx_to_rid = self.kegg.create_reaction_dicts()
        self.kegg.cid_to_idx, self.kegg.idx_to_cid = self.kegg.create_compound_dicts()
        self.kegg.S = self.kegg.create_S_from_irreversible_network()
        x0 = self.kegg.initialize_metabolite_vector(self.seedSet)
        R = (self.kegg.S < 0)*1
        P = (self.kegg.S > 0)*1
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        x0 = csr_matrix(x0)
        x0 = x0.transpose()

        self.X,self.Y = ne.netExp_trace(R,P,x0,b)
    
    def test_iteration_dict(self):
        compound_iteration_dict = self.kegg.create_iteration_dict(self.X, self.kegg.idx_to_cid)
        reaction_iteration_dict = self.kegg.create_iteration_dict(self.Y, self.kegg.idx_to_rid)

        ## Look through iteration 0 compounds and make sure its exactly equal to seed set
        self.assertEqual(set([k for k,v in compound_iteration_dict.items() if v==0]), set(self.seedSet))

        ## Check that matrix has exactlly two more rows than max(idx_iter.values())
        self.assertEqual(len(self.X)-2, max(compound_iteration_dict.values()) )
        self.assertEqual(len(self.Y)-1, max(reaction_iteration_dict.values()) )
        ## Check that all iterations are present in at least 1 compound in the iteration dict
        self.assertEqual(len(self.X)-1, len(set(compound_iteration_dict.values())) )
        self.assertEqual(len(self.Y)-1, len(set(reaction_iteration_dict.values())) )


class TestLoadTupleNetwork(unittest.TestCase):

    def test_tuple_len_throws(self):

        with self.assertRaises(ValueError):
            ne._load_tuple_network([([1])])

        with self.assertRaises(ValueError):
            ne._load_tuple_network([([1],[2],[3])])

        rxns = [
             (["A","B"],["C"]),
             (["C","D"],["E","F"]),
             (["E","F"],["G"]),
             (["G","H"],["I"]),
             (["A","J"],["I"])]

        expected_df = pd.DataFrame(
            [{'rn': 0, 'cid': 'A', 's': -1},
            {'rn': 0, 'cid': 'B', 's': -1},
            {'rn': 0, 'cid': 'C', 's': 1},
            {'rn': 1, 'cid': 'C', 's': -1},
            {'rn': 1, 'cid': 'D', 's': -1},
            {'rn': 1, 'cid': 'E', 's': 1},
            {'rn': 1, 'cid': 'F', 's': 1},
            {'rn': 2, 'cid': 'E', 's': -1},
            {'rn': 2, 'cid': 'F', 's': -1},
            {'rn': 2, 'cid': 'G', 's': 1},
            {'rn': 3, 'cid': 'G', 's': -1},
            {'rn': 3, 'cid': 'H', 's': -1},
            {'rn': 3, 'cid': 'I', 's': 1},
            {'rn': 4, 'cid': 'A', 's': -1},
            {'rn': 4, 'cid': 'J', 's': -1},
            {'rn': 4, 'cid': 'I', 's': 1}])

        self.assertTrue(dfs_have_equal_content(ne._load_tuple_network(rxns),expected_df))


class TestGlobalMetabolicNetworkExpand(unittest.TestCase):

    def setUp(self):
        self.toy = ne.GlobalMetabolicNetwork("dev")
        rnxs = [
            (["A","B"],["C"]),
            (["C","D"],["E","F"]),
            (["E","F"],["G"]),
            (["G","H"],["I"]),
            (["A","J"],["I"])]

        self.toy.network = ne._load_tuple_network(rnxs)
        self.toy.convertToIrreversible()

    def test_expansion_1(self):

        compounds, reactions = self.toy.expand(["A","B","D","H"],"trace")

        expected_compounds = {'H': 0,
                            'B': 0,
                            'D': 0,
                            'A': 0,
                            'C': 1,
                            'F': 2,
                            'E': 2,
                            'G': 3,
                            'I': 4,
                            'J': 5}

        expected_reactions = {(0, 'forward'): 1,
                            (0, 'reverse'): 2,
                            (1, 'forward'): 2,
                            (2, 'forward'): 3,
                            (1, 'reverse'): 3,
                            (2, 'reverse'): 4,
                            (3, 'forward'): 4,
                            (4, 'reverse'): 5,
                            (3, 'reverse'): 5,
                            (4, 'forward'): 6}

        self.assertEqual(compounds,expected_compounds)
        self.assertEqual(reactions,expected_reactions)