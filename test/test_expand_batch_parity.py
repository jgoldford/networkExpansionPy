"""
Parity tests for the consolidated Rust API:
  - expand_batch (unified, replaces expand/expand_masked/expand_masked_batch)
  - expand_trace_batch (new Rust trace)
  - contract_batch (unified, replaces contract/contract_batch)
"""

import unittest
import numpy as np
from random import sample, seed as random_seed

import networkExpansionPy.lib as ne

try:
    import netexprs
except ImportError:
    netexprs = None


# ── Helpers ──────────────────────────────────────────────────────────────

def make_toy():
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
    return toy


# ── expand_batch: low-level ──────────────────────────────────────────────

@unittest.skipIf(netexprs is None, "netexprs not installed")
class TestExpandBatchLowLevel(unittest.TestCase):
    """Direct netexprs.expand_batch calls — compare 1-row vs N-row consistency."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        cls.kegg._ensure_rust_ready()
        cls.ra = cls.kegg._rust_arrays
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.all_seeds = list(set(seeds_df["ID"]))

    def _batch_call(self, x_batch, masks=None):
        ra = self.ra
        return netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], ra["n_reactions"],
            ra["p_data"], ra["p_indices"], ra["p_indptr"], ra["n_compounds"],
            x_batch, ra["b"], masks,
        )

    def test_single_row_no_mask(self):
        x0 = self.kegg.initialize_metabolite_vector(self.all_seeds).astype(np.uint8)
        x_out, y_out = self._batch_call(x0.reshape(1, -1))
        self.assertEqual(x_out.shape[0], 1)
        self.assertTrue(x_out[0].sum() > 0)

    def test_multi_row_no_mask(self):
        random_seed(42)
        seed_sets = [self.all_seeds, self.all_seeds[:10], sample(self.all_seeds, 20)]
        n = len(seed_sets)
        x_batch = np.zeros((n, self.ra["n_compounds"]), dtype=np.uint8)
        for i, s in enumerate(seed_sets):
            x_batch[i] = self.kegg.initialize_metabolite_vector(s).astype(np.uint8)

        x_out, y_out = self._batch_call(x_batch)

        # Each row should match a single-row call
        for i in range(n):
            x_single, y_single = self._batch_call(x_batch[i:i+1])
            np.testing.assert_array_equal(x_out[i], x_single[0], err_msg=f"x mismatch row {i}")
            np.testing.assert_array_equal(y_out[i], y_single[0], err_msg=f"y mismatch row {i}")

    def test_with_mask(self):
        """Single seed + mask should match old expand_masked behavior."""
        x0 = self.kegg.initialize_metabolite_vector(self.all_seeds).astype(np.uint8)
        mask = np.ones(self.ra["n_reactions"], dtype=np.uint8)
        # Zero out first 100 reactions
        mask[:100] = 0

        x_out, y_out = self._batch_call(x0.reshape(1, -1), mask.reshape(1, -1))
        # Masked-out reactions should not fire
        np.testing.assert_array_equal(y_out[0, :100], 0)

    def test_multi_row_with_masks(self):
        """N seeds × N masks — each row gets its own mask."""
        random_seed(99)
        n = 5
        x_batch = np.zeros((n, self.ra["n_compounds"]), dtype=np.uint8)
        masks = np.ones((n, self.ra["n_reactions"]), dtype=np.uint8)
        for i in range(n):
            seeds = sample(self.all_seeds, 10 + i * 5)
            x_batch[i] = self.kegg.initialize_metabolite_vector(seeds).astype(np.uint8)
            masks[i, :i*50] = 0  # increasingly aggressive masking

        x_out_batch, y_out_batch = self._batch_call(x_batch, masks)

        for i in range(n):
            x_single, y_single = self._batch_call(x_batch[i:i+1], masks[i:i+1])
            np.testing.assert_array_equal(x_out_batch[i], x_single[0])
            np.testing.assert_array_equal(y_out_batch[i], y_single[0])

    def test_empty_seeds(self):
        """All-zero seed row should produce empty expansion."""
        x_batch = np.zeros((1, self.ra["n_compounds"]), dtype=np.uint8)
        x_out, y_out = self._batch_call(x_batch)
        self.assertEqual(x_out[0].sum(), 0)
        self.assertEqual(y_out[0].sum(), 0)


# ── expand_batch: high-level ─────────────────────────────────────────────

class TestExpandHighLevel(unittest.TestCase):
    """_expand_rust and run_expansions_batch vs Python equivalents."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.all_seeds = list(set(seeds_df["ID"]))

    def test_single_expand_parity(self):
        cpds_r, rxns_r = self.kegg._expand_rust(self.all_seeds)
        cpds_p, rxns_p = self.kegg._expand_python(self.all_seeds, "naive")
        self.assertEqual(set(cpds_r), set(cpds_p))
        self.assertEqual(set(rxns_r), set(rxns_p))

    def test_single_expand_with_mask(self):
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())
        to_exclude = sample(all_rxns, len(all_rxns) // 10)

        cpds_r, rxns_r = self.kegg._expand_rust(self.all_seeds, excluded_reactions=to_exclude)
        cpds_p, rxns_p = self.kegg._expand_python(self.all_seeds, "naive", excluded_reactions=to_exclude)
        self.assertEqual(set(cpds_r), set(cpds_p))
        self.assertEqual(set(rxns_r), set(rxns_p))

    def test_run_expansions_batch_parity(self):
        random_seed(42)
        seed_sets = [self.all_seeds, self.all_seeds[:10], sample(self.all_seeds, 20)]

        cpds_r, rxns_r = self.kegg.run_expansions_batch(seed_sets, algorithm="naive")
        cpds_p, rxns_p = self.kegg._run_expansions_python(seed_sets, "naive")

        for i in range(len(seed_sets)):
            self.assertEqual(set(cpds_r[i]), set(cpds_p[i]), f"cpds mismatch set {i}")
            self.assertEqual(set(rxns_r[i]), set(rxns_p[i]), f"rxns mismatch set {i}")

    def test_run_expansions_reactionMasks_parity(self):
        random_seed(123)
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())
        network_rxns = list(set(r[0] for r in all_rxns))

        masked_sets = []
        for _ in range(10):
            n_remove = np.random.randint(0, len(network_rxns) // 10)
            base = sample(network_rxns, n_remove)
            masked_sets.append([r for r in all_rxns if r[0] in base])

        cpds_r, rxns_r = self.kegg.run_expansions_batch(self.all_seeds, masked_sets)
        cpds_p, rxns_p = self.kegg._run_expansions_reactionMasks_python(self.all_seeds, masked_sets)

        for i in range(len(masked_sets)):
            self.assertEqual(set(cpds_r[i]), set(cpds_p[i]), f"cpds mismatch mask {i}")
            self.assertEqual(set(rxns_r[i]), set(rxns_p[i]), f"rxns mismatch mask {i}")


# ── expand_trace_batch ───────────────────────────────────────────────────

class TestExpandTraceParity(unittest.TestCase):
    """Rust trace vs Python trace."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.all_seeds = list(set(seeds_df["ID"]))

    def test_single_trace_parity(self):
        cpds_r, rxns_r = self.kegg._expand_trace_rust(self.all_seeds)
        cpds_p, rxns_p = self.kegg._expand_python(self.all_seeds, "trace")
        self.assertEqual(cpds_r, cpds_p)
        self.assertEqual(rxns_r, rxns_p)

    def test_single_trace_small_seed(self):
        cpds_r, rxns_r = self.kegg._expand_trace_rust(["C00001"])
        cpds_p, rxns_p = self.kegg._expand_python(["C00001"], "trace")
        self.assertEqual(cpds_r, cpds_p)
        self.assertEqual(rxns_r, rxns_p)

    def test_batch_trace_parity(self):
        random_seed(42)
        seed_sets = [self.all_seeds, self.all_seeds[:10], sample(self.all_seeds, 20)]

        cpds_r, rxns_r = self.kegg.run_expansions_batch(seed_sets, algorithm="trace")
        cpds_p, rxns_p = self.kegg._run_expansions_python(seed_sets, "trace")

        for i in range(len(seed_sets)):
            self.assertEqual(cpds_r[i], cpds_p[i], f"compound trace mismatch set {i}")
            self.assertEqual(rxns_r[i], rxns_p[i], f"reaction trace mismatch set {i}")

    def test_trace_via_public_api(self):
        """expand() with algorithm='trace' should dispatch to Rust and match Python."""
        cpds_api, rxns_api = self.kegg.expand(self.all_seeds, algorithm="trace")
        cpds_py, rxns_py = self.kegg._expand_python(self.all_seeds, "trace")
        self.assertEqual(cpds_api, cpds_py)
        self.assertEqual(rxns_api, rxns_py)


class TestExpandTraceParityToy(unittest.TestCase):
    """Trace tests on toy network — easy to debug."""

    @classmethod
    def setUpClass(cls):
        cls.toy = make_toy()

    def test_toy_trace_full_seeds(self):
        cpds_r, rxns_r = self.toy._expand_trace_rust(["A", "B", "D", "H"])
        cpds_p, rxns_p = self.toy._expand_python(["A", "B", "D", "H"], "trace")
        self.assertEqual(cpds_r, cpds_p)
        self.assertEqual(rxns_r, rxns_p)

    def test_toy_trace_single_seed(self):
        cpds_r, rxns_r = self.toy._expand_trace_rust(["I"])
        cpds_p, rxns_p = self.toy._expand_python(["I"], "trace")
        self.assertEqual(cpds_r, cpds_p)
        self.assertEqual(rxns_r, rxns_p)

    def test_toy_trace_no_expansion(self):
        """Seeds that can't fire anything — trace should only have iteration 0."""
        cpds_r, rxns_r = self.toy._expand_trace_rust(["A", "D", "E", "H"])
        cpds_p, rxns_p = self.toy._expand_python(["A", "D", "E", "H"], "trace")
        self.assertEqual(cpds_r, cpds_p)
        self.assertEqual(rxns_r, rxns_p)
        # All compounds at iteration 0, no reactions
        self.assertTrue(all(v == 0 for v in cpds_r.values()))
        self.assertEqual(len(rxns_r), 0)


# ── contract_batch: high-level ───────────────────────────────────────────

class TestContractHighLevel(unittest.TestCase):
    """Consolidated contract_batch vs Python contraction."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.all_seeds = list(set(seeds_df["ID"]))
        cls.cpds_scope, cls.rxns_scope = cls.kegg._expand_python(cls.all_seeds, "naive")

    def test_single_contraction_parity(self):
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())
        extinct = sample(all_rxns, len(all_rxns) // 20)

        cpds_r, rxns_r = self.kegg._contract_rust(self.rxns_scope, self.cpds_scope, extinct)
        cpds_p, rxns_p = self.kegg._contract_python(self.rxns_scope, self.cpds_scope, extinct)
        self.assertEqual(set(cpds_r), set(cpds_p))
        self.assertEqual(set(rxns_r), set(rxns_p))

    def test_batch_contraction_parity(self):
        random_seed(42)
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())
        network_rxns = list(set(r[0] for r in all_rxns))

        extinct_sets = []
        for _ in range(10):
            n = np.random.randint(1, len(network_rxns) // 10)
            base = sample(network_rxns, n)
            extinct_sets.append([r for r in all_rxns if r[0] in base])

        cpds_r, rxns_r = self.kegg._run_contractions_rust(
            self.rxns_scope, self.cpds_scope, extinct_sets)
        cpds_p, rxns_p = self.kegg._run_contractions_python(
            self.rxns_scope, self.cpds_scope, extinct_sets)

        for i in range(len(extinct_sets)):
            self.assertEqual(set(cpds_r[i]), set(cpds_p[i]), f"cpds mismatch batch {i}")
            self.assertEqual(set(rxns_r[i]), set(rxns_p[i]), f"rxns mismatch batch {i}")


# ── Dispatch tests ───────────────────────────────────────────────────────

class TestDispatch(unittest.TestCase):
    """Verify public methods route to Rust when available."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.all_seeds = list(set(seeds_df["ID"]))

    @unittest.skipIf(netexprs is None, "netexprs not installed")
    def test_expand_naive_dispatches_rust(self):
        from unittest.mock import patch
        with patch.object(self.kegg, '_expand_rust', wraps=self.kegg._expand_rust) as m:
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.expand(self.all_seeds[:5], algorithm="naive")
            m.assert_called_once()

    @unittest.skipIf(netexprs is None, "netexprs not installed")
    def test_expand_trace_dispatches_rust(self):
        from unittest.mock import patch
        with patch.object(self.kegg, '_expand_trace_rust', wraps=self.kegg._expand_trace_rust) as m:
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.expand(self.all_seeds[:5], algorithm="trace")
            m.assert_called_once()

    def test_expand_cr_stays_python(self):
        """'cr' and 'step' algorithms should always use Python."""
        from unittest.mock import patch
        with patch.object(self.kegg, '_expand_python', wraps=self.kegg._expand_python) as m:
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.expand(self.all_seeds[:5], algorithm="cr")
            m.assert_called_once()


if __name__ == "__main__":
    unittest.main()