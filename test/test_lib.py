import unittest
import networkExpansionPy.lib as ne

class TestGlobalMetabolicNetworkInit(unittest.TestCase):

    def test_kegg_og_no_args(self):
        kegg_og = ne.GlobalMetabolicNetwork()
        kegg_og.pruneUnbalancedReactions()

    def test_kegg_og(self):
        kegg_og = ne.GlobalMetabolicNetwork("KEGG_OG")

    def test_ecg(self):
        ecg = ne.GlobalMetabolicNetwork("ecg")

    def test_KEGG(self):
        kegg = ne.GlobalMetabolicNetwork("KEGG")