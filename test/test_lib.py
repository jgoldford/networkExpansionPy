import unittest
import networkExpansionPy.lib as ne
import pandas as pd

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