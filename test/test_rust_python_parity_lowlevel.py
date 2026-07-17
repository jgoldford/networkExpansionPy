import unittest
import numpy as np
from scipy.sparse import csr_matrix
import networkExpansionPy.lib as ne
import netexprs
import pandas as pd

def rpxb(m, seedSet):
    """Extract R, P, x0, b matrices from a GlobalMetabolicNetwork."""
    m._ensure_dicts()
    x0 = m.initialize_metabolite_vector(seedSet)
    R, P = m.create_RP_from_irreversible_network()
    b = sum(R)

    R = csr_matrix(R)
    P = csr_matrix(P)
    b = csr_matrix(b).transpose()
    x0 = csr_matrix(x0).transpose()

    return R, P, x0, b


def netExp_rs(R, P, x, b):
    """Rust wrapper matching netExp signature — 1-row batch, no mask."""
    R_T = R.T.tocsr()
    P = P.tocsr()

    x_init = np.asarray(x.toarray()).ravel().astype(np.uint8).reshape(1, -1)

    x_out, y_out = netexprs.expand_batch(
        R_T.data.astype(np.float64),
        R_T.indices.astype(np.int32),
        R_T.indptr.astype(np.int32),
        R_T.shape[0],
        P.data.astype(np.float64),
        P.indices.astype(np.int32),
        P.indptr.astype(np.int32),
        P.shape[0],
        x_init,
        np.asarray(b.toarray()).ravel().astype(np.float64),
        None,
    )

    return x_out[0], y_out[0]

def netExp_masked_rs(R, P, x, b, mask):
    """Rust wrapper for masked expansion — 1-row batch with 1-row mask."""
    R_T = R.T.tocsr()
    P = P.tocsr()

    x_init = np.asarray(x.toarray()).ravel().astype(np.uint8).reshape(1, -1)
    mask_2d = mask.astype(np.uint8).reshape(1, -1)

    x_out, y_out = netexprs.expand_batch(
        R_T.data.astype(np.float64),
        R_T.indices.astype(np.int32),
        R_T.indptr.astype(np.int32),
        R_T.shape[0],
        P.data.astype(np.float64),
        P.indices.astype(np.int32),
        P.indptr.astype(np.int32),
        P.shape[0],
        x_init,
        np.asarray(b.toarray()).ravel().astype(np.float64),
        mask_2d,
    )

    return x_out[0], y_out[0]

def netExp_masked_py(R, P, x, b, mask):
    """Python masked expansion using diagonal matrix approach."""
    reaction_mask = csr_matrix(np.diag(mask))
    Rstar = R * reaction_mask
    Pstar = P * reaction_mask
    return ne.netExp(Rstar, Pstar, x, b)

def netExp_masked_batch_rs(R, P, x, b, masks):
    """Rust wrapper for batch masked expansion — tiled seed, N masks."""
    R_T = R.T.tocsr()
    P = P.tocsr()

    x_1d = np.asarray(x.toarray()).ravel().astype(np.uint8)
    x_init = np.tile(x_1d, (masks.shape[0], 1))

    x_out, y_out = netexprs.expand_batch(
        R_T.data.astype(np.float64),
        R_T.indices.astype(np.int32),
        R_T.indptr.astype(np.int32),
        R_T.shape[0],
        P.data.astype(np.float64),
        P.indices.astype(np.int32),
        P.indptr.astype(np.int32),
        P.shape[0],
        x_init,
        np.asarray(b.toarray()).ravel().astype(np.float64),
        masks.astype(np.uint8),
    )

    return x_out, y_out

def netContract_rs(R, P, x_active, y_active, y_extinct):
    """Rust wrapper for network contraction — 1-row batch."""
    R_T = R.T.tocsr()
    P = P.tocsr()

    x_out, y_out = netexprs.contract_batch(
        R_T.data.astype(np.float64),
        R_T.indices.astype(np.int32),
        R_T.indptr.astype(np.int32),
        R_T.shape[0],
        P.data.astype(np.float64),
        P.indices.astype(np.int32),
        P.indptr.astype(np.int32),
        P.shape[0],
        np.asarray(x_active.toarray()).ravel().astype(np.uint8).reshape(1, -1),
        np.asarray(y_active.toarray()).ravel().astype(np.uint8).reshape(1, -1),
        np.asarray(y_extinct.toarray()).ravel().astype(np.uint8).reshape(1, -1),
        None,
    )
    return x_out[0], y_out[0]

def netContract_batch_rs(R, P, x_active, y_active, y_extinct_batch):
    """Rust wrapper for batch network contraction — tiled scope, N extinction sets."""
    R_T = R.T.tocsr()
    P = P.tocsr()

    n_batches = y_extinct_batch.shape[0]
    x_1d = np.asarray(x_active.toarray()).ravel().astype(np.uint8)
    y_1d = np.asarray(y_active.toarray()).ravel().astype(np.uint8)

    x_out, y_out = netexprs.contract_batch(
        R_T.data.astype(np.float64),
        R_T.indices.astype(np.int32),
        R_T.indptr.astype(np.int32),
        R_T.shape[0],
        P.data.astype(np.float64),
        P.indices.astype(np.int32),
        P.indptr.astype(np.int32),
        P.shape[0],
        np.tile(x_1d, (n_batches, 1)),
        np.tile(y_1d, (n_batches, 1)),
        y_extinct_batch.astype(np.uint8),
        None,
    )
    return x_out, y_out

class TestExpandParity(unittest.TestCase):

    def test_toy_network(self):
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        x_py, y_py = ne.netExp(R, P, x0, b)
        x_rs, y_rs = netExp_rs(R, P, x0, b)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_toy_network_different_seed(self):
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["I"])

        x_py, y_py = ne.netExp(R, P, x0, b)
        x_rs, y_rs = netExp_rs(R, P, x0, b)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_toy_network_no_expansion(self):
        """Seed set that cannot fire any reactions."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "D", "E", "H"])

        x_py, y_py = ne.netExp(R, P, x0, b)
        x_rs, y_rs = netExp_rs(R, P, x0, b)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())
class TestExpandMaskedParity(unittest.TestCase):

    def test_toy_network_masked(self):
        """Mask out one reaction."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        # Mask out reaction 0 (forward and reverse)
        n_reactions = R.shape[1]
        mask = np.ones(n_reactions, dtype=np.uint8)
        mask[0] = 0
        mask[1] = 0

        x_py, y_py = netExp_masked_py(R, P, x0, b, mask)
        x_rs, y_rs = netExp_masked_rs(R, P, x0, b, mask)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_toy_network_half_masked(self):
        """Mask out half the reactions."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        n_reactions = R.shape[1]
        mask = np.ones(n_reactions, dtype=np.uint8)
        mask[: n_reactions // 2] = 0

        x_py, y_py = netExp_masked_py(R, P, x0, b, mask)
        x_rs, y_rs = netExp_masked_rs(R, P, x0, b, mask)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())
class TestExpandLargeNetwork(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Load KEGG network once for all tests in this class."""
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_network(self):
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        x_py, y_py = ne.netExp(R, P, x0, b)
        x_rs, y_rs = netExp_rs(R, P, x0, b)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_kegg_network_minimal_seed(self):
        """Single compound seed."""
        seedSet = ["C00001"]

        R, P, x0, b = rpxb(self.kegg, seedSet)

        x_py, y_py = ne.netExp(R, P, x0, b)
        x_rs, y_rs = netExp_rs(R, P, x0, b)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())


class TestExpandMaskedLargeNetwork(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Load KEGG network once for all tests in this class."""
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_masked_random(self):
        """Randomly mask 10% of reactions."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        n_reactions = R.shape[1]
        np.random.seed(42)
        mask = (np.random.random(n_reactions) > 0.1).astype(np.uint8)

        x_py, y_py = netExp_masked_py(R, P, x0, b, mask)
        x_rs, y_rs = netExp_masked_rs(R, P, x0, b, mask)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_kegg_masked_half(self):
        """Mask out first half of reactions."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        n_reactions = R.shape[1]
        mask = np.ones(n_reactions, dtype=np.uint8)
        mask[: n_reactions // 2] = 0

        x_py, y_py = netExp_masked_py(R, P, x0, b, mask)
        x_rs, y_rs = netExp_masked_rs(R, P, x0, b, mask)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

class TestExpandMaskedBatchParity(unittest.TestCase):

    def test_toy_batch(self):
        """Run multiple masked expansions in batch."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])
        n_reactions = R.shape[1]

        # Create 5 different masks
        np.random.seed(123)
        masks = (np.random.random((5, n_reactions)) > 0.3).astype(np.uint8)

        x_batch, y_batch = netExp_masked_batch_rs(R, P, x0, b, masks)

        # Compare each row against individual Python runs
        for i in range(masks.shape[0]):
            x_py, y_py = netExp_masked_py(R, P, x0, b, masks[i])
            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

    @classmethod
    def setUpClass(cls):
        """Load real metabolic network once for all tests in this class."""
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_batch(self):
        """Batch masked expansion on real network."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)
        n_reactions = R.shape[1]

        # Create 20 different masks
        np.random.seed(456)
        masks = (np.random.random((20, n_reactions)) > 0.1).astype(np.uint8)

        x_batch, y_batch = netExp_masked_batch_rs(R, P, x0, b, masks)

        # Compare each row against individual Python runs
        for i in range(masks.shape[0]):
            x_py, y_py = netExp_masked_py(R, P, x0, b, masks[i])
            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

    def test_kegg_batch_large(self):
        """Batch masked expansion with many masks."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)
        n_reactions = R.shape[1]

        # Create 100 different masks
        np.random.seed(789)
        masks = (np.random.random((50, n_reactions)) > 0.1).astype(np.uint8)

        x_batch, y_batch = netExp_masked_batch_rs(R, P, x0, b, masks)

        # Compare each row against individual Python runs
        for i in range(masks.shape[0]):
            x_py, y_py = netExp_masked_py(R, P, x0, b, masks[i])
            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

class TestContractParity(unittest.TestCase):

    def test_toy_contract(self):
        """Contract after expansion, removing one reaction."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        # First expand to get full scope
        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        # Mark some reactions as extinct
        n_reactions = R.shape[1]
        y_extinct = csr_matrix(np.zeros(n_reactions)).T
        y_extinct[0, 0] = 1  # Kill first reaction

        # Python contraction
        Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
        x_py = Xac_py[-1]
        y_py = Yac_py[-1]

        # Rust contraction
        x_rs, y_rs = netContract_rs(R, P, x_expanded, y_expanded, y_extinct)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_toy_contract_multiple_extinct(self):
        """Contract with multiple extinct reactions."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        y_extinct = csr_matrix(np.zeros(n_reactions)).T
        y_extinct[0, 0] = 1
        y_extinct[2, 0] = 1
        y_extinct[4, 0] = 1

        Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
        x_py = Xac_py[-1]
        y_py = Yac_py[-1]

        x_rs, y_rs = netContract_rs(R, P, x_expanded, y_expanded, y_extinct)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_toy_contract_no_extinction(self):
        """Contract with no extinct reactions (should be no-op)."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        y_extinct = csr_matrix(np.zeros(n_reactions)).T

        Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
        x_py = Xac_py[-1]
        y_py = Yac_py[-1]

        x_rs, y_rs = netContract_rs(R, P, x_expanded, y_expanded, y_extinct)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())


class TestContractLargeNetwork(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_contract(self):
        """Contract on real network."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(42)
        y_extinct_arr = (np.random.random(n_reactions) < 0.05).astype(np.int32)
        y_extinct = csr_matrix(y_extinct_arr).T

        Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
        x_py = Xac_py[-1]
        y_py = Yac_py[-1]

        x_rs, y_rs = netContract_rs(R, P, x_expanded, y_expanded, y_extinct)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

    def test_kegg_contract_heavy_extinction(self):
        """Contract with many extinct reactions."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(123)
        y_extinct_arr = (np.random.random(n_reactions) < 0.3).astype(np.int32)
        y_extinct = csr_matrix(y_extinct_arr).T

        Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
        x_py = Xac_py[-1]
        y_py = Yac_py[-1]

        x_rs, y_rs = netContract_rs(R, P, x_expanded, y_expanded, y_extinct)

        np.testing.assert_array_equal(x_rs, x_py.toarray().ravel())
        np.testing.assert_array_equal(y_rs, y_py.toarray().ravel())

class TestContractBatchParity(unittest.TestCase):

    def test_toy_contract_batch(self):
        """Batch contraction on toy network."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(42)
        n_batches = 5
        y_extinct_batch = (np.random.random((n_batches, n_reactions)) < 0.3).astype(np.uint8)

        x_batch, y_batch = netContract_batch_rs(R, P, x_expanded, y_expanded, y_extinct_batch)

        for i in range(n_batches):
            y_extinct = csr_matrix(y_extinct_batch[i]).T
            Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
            x_py = Xac_py[-1]
            y_py = Yac_py[-1]

            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

    def test_toy_contract_batch_varying_extinction(self):
        """Batch contraction with varying extinction rates."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
            (["G", "H"], ["I"]),
            (["A", "J"], ["I"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        R, P, x0, b = rpxb(toy, ["A", "B", "D", "H"])

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(123)

        # Different extinction rates per batch
        y_extinct_batch = np.vstack([
            (np.random.random(n_reactions) < 0.1).astype(np.uint8),
            (np.random.random(n_reactions) < 0.3).astype(np.uint8),
            (np.random.random(n_reactions) < 0.5).astype(np.uint8),
            np.zeros(n_reactions, dtype=np.uint8),
            np.ones(n_reactions, dtype=np.uint8),
        ])

        x_batch, y_batch = netContract_batch_rs(R, P, x_expanded, y_expanded, y_extinct_batch)

        for i in range(y_extinct_batch.shape[0]):
            y_extinct = csr_matrix(y_extinct_batch[i]).T
            Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
            x_py = Xac_py[-1]
            y_py = Yac_py[-1]
            
            print(f"Batch {i}: extinct_sum={y_extinct_batch[i].sum()}")
            print(f"  x_rs: {x_batch[i].sum()}, x_py: {x_py.toarray().sum()}")
            print(f"  y_rs: {y_batch[i].sum()}, y_py: {y_py.toarray().sum()}")

            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())


class TestContractBatchLargeNetwork(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_contract_batch(self):
        """Batch contraction on real network."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(456)
        n_batches = 20
        y_extinct_batch = (np.random.random((n_batches, n_reactions)) < 0.1).astype(np.uint8)

        x_batch, y_batch = netContract_batch_rs(R, P, x_expanded, y_expanded, y_extinct_batch)

        for i in range(n_batches):
            y_extinct = csr_matrix(y_extinct_batch[i]).T
            Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
            x_py = Xac_py[-1]
            y_py = Yac_py[-1]

            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

    def test_kegg_contract_batch_large(self):
        """Batch contraction with many extinction sets."""
        R, P, x0, b = rpxb(self.kegg, self.seedSet)

        x_expanded, y_expanded = ne.netExp(R, P, x0, b)

        n_reactions = R.shape[1]
        np.random.seed(789)
        n_batches = 100
        y_extinct_batch = (np.random.random((n_batches, n_reactions)) < 0.1).astype(np.uint8)

        x_batch, y_batch = netContract_batch_rs(R, P, x_expanded, y_expanded, y_extinct_batch)

        for i in range(n_batches):
            y_extinct = csr_matrix(y_extinct_batch[i]).T
            Xac_py, Yac_py = ne.netContract(R, P, b, x_expanded, y_expanded, y_extinct)
            x_py = Xac_py[-1]
            y_py = Yac_py[-1]

            np.testing.assert_array_equal(x_batch[i], x_py.toarray().ravel())
            np.testing.assert_array_equal(y_batch[i], y_py.toarray().ravel())

if __name__ == "__main__":
    unittest.main()