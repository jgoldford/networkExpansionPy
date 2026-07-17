"""
Tests for batch initialization methods and contract_and_reexpand().

Tests cover:
1. initialize_reaction_matrix — batch version of initialize_reaction_vector
2. initialize_metabolite_matrix — batch version of initialize_metabolite_vector
3. contract_and_reexpand — fused contraction + re-expansion without ID marshalling
"""

import unittest
import numpy as np
from random import sample, seed as random_seed
import networkExpansionPy.lib as ne


class TestInitializeReactionMatrix(unittest.TestCase):
    """Tests for batch reaction vector initialization."""

    def setUp(self):
        self.toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
        ]
        self.toy.network = ne._load_tuple_network(rxns)
        self.toy.convertToIrreversible()
        self.toy._ensure_dicts()

    def test_matches_individual_calls(self):
        """Matrix rows must match individual initialize_reaction_vector calls."""
        all_rxns = list(self.toy.rid_to_idx.keys())
        sets = [
            all_rxns[:2],
            all_rxns[2:4],
            all_rxns,
            [],
        ]

        mat = self.toy.initialize_reaction_matrix(sets)

        self.assertEqual(mat.shape, (4, len(self.toy.rid_to_idx)))
        self.assertEqual(mat.dtype, np.uint8)

        for i, s in enumerate(sets):
            expected = self.toy.initialize_reaction_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)

    def test_string_ids_expand_both_directions(self):
        """Plain rn ID strings should expand to both forward and reverse."""
        # Use plain string IDs (not tuples)
        string_sets = [["0"], ["1", "2"]]
        mat = self.toy.initialize_reaction_matrix(string_sets)

        for i, s in enumerate(string_sets):
            expected = self.toy.initialize_reaction_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)

    def test_empty_input(self):
        """Empty list of sets returns empty matrix."""
        mat = self.toy.initialize_reaction_matrix([])
        self.assertEqual(mat.shape, (0, len(self.toy.rid_to_idx)))


class TestInitializeReactionMatrixKEGG(unittest.TestCase):
    """Tests on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.kegg._ensure_dicts()
        cls.all_rxns = list(cls.kegg.rid_to_idx.keys())

    def test_kegg_matches_individual_calls(self):
        """Matrix rows match individual calls on real network."""
        random_seed(42)
        sets = [sample(self.all_rxns, 100) for _ in range(10)]

        mat = self.kegg.initialize_reaction_matrix(sets)

        for i, s in enumerate(sets):
            expected = self.kegg.initialize_reaction_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)

    def test_kegg_string_ids(self):
        """Plain rn strings on KEGG network."""
        # Get unique base rn IDs (strings, no direction)
        base_rns = list(set(r[0] for r in self.all_rxns))[:50]
        sets = [base_rns[:10], base_rns[10:30], base_rns]

        mat = self.kegg.initialize_reaction_matrix(sets)

        for i, s in enumerate(sets):
            expected = self.kegg.initialize_reaction_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)


class TestInitializeMetaboliteMatrix(unittest.TestCase):
    """Tests for batch metabolite vector initialization."""

    def setUp(self):
        self.toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
        ]
        self.toy.network = ne._load_tuple_network(rxns)
        self.toy.convertToIrreversible()
        self.toy._ensure_dicts()

    def test_matches_individual_calls(self):
        """Matrix rows must match individual initialize_metabolite_vector calls."""
        sets = [
            ["A", "B"],
            ["C", "D", "E"],
            ["A", "B", "C", "D", "E", "F", "G"],
            [],
        ]

        mat = self.toy.initialize_metabolite_matrix(sets)

        self.assertEqual(mat.shape, (4, len(self.toy.cid_to_idx)))
        self.assertEqual(mat.dtype, np.uint8)

        for i, s in enumerate(sets):
            expected = self.toy.initialize_metabolite_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)

    def test_unknown_ids_ignored(self):
        """Unknown compound IDs are silently ignored."""
        mat = self.toy.initialize_metabolite_matrix([["A", "UNKNOWN", "B"]])
        expected = self.toy.initialize_metabolite_vector(["A", "B"]).astype(np.uint8)
        np.testing.assert_array_equal(mat[0], expected)


class TestInitializeMetaboliteMatrixKEGG(unittest.TestCase):
    """Tests on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.kegg._ensure_dicts()
        cls.seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_matches_individual_calls(self):
        """Matrix rows match individual calls on real network."""
        random_seed(42)
        all_cpds = list(self.kegg.cid_to_idx.keys())
        sets = [sample(all_cpds, 50) for _ in range(10)]

        mat = self.kegg.initialize_metabolite_matrix(sets)

        for i, s in enumerate(sets):
            expected = self.kegg.initialize_metabolite_vector(s).astype(np.uint8)
            np.testing.assert_array_equal(mat[i], expected)


class TestContractAndReexpand(unittest.TestCase):
    """Tests for the fused contract + re-expand method."""

    def test_toy_matches_separate_calls(self):
        """contract_and_reexpand must match separate contraction + expansion."""
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
        toy._ensure_dicts()

        seeds = ["A", "B", "D", "H"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        all_rxns = list(toy.rid_to_idx.keys())
        extinct_sets = [all_rxns[:2], all_rxns[2:5], all_rxns[::3]]

        result = toy.contract_and_reexpand(seeds, rxns_scope, cpds_scope, extinct_sets)

        # Compare against separate calls
        cpds_c, rxns_c = toy._run_contractions_rust(rxns_scope, cpds_scope, extinct_sets, seeds)
        for i in range(len(extinct_sets)):
            # Check contraction output
            expected_x = toy.initialize_metabolite_vector(cpds_c[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['x_contracted'][i], expected_x)

            expected_y = toy.initialize_reaction_vector(rxns_c[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['y_contracted'][i], expected_y)

        # Check re-expansion output against run_expansions_batch
        cpds_r, rxns_r = toy.run_expansions_batch(cpds_c, extinct_sets)
        for i in range(len(extinct_sets)):
            expected_x = toy.initialize_metabolite_vector(cpds_r[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['x_reexpanded'][i], expected_x)

            expected_y = toy.initialize_reaction_vector(rxns_r[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['y_reexpanded'][i], expected_y)

    def test_counts_from_raw_arrays(self):
        """Demonstrate extracting counts directly from raw arrays."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
            (["E", "F"], ["G"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B", "D"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        all_rxns = list(toy.rid_to_idx.keys())
        extinct_sets = [all_rxns[:2], all_rxns[2:4]]

        result = toy.contract_and_reexpand(seeds, rxns_scope, cpds_scope, extinct_sets)

        # Counts directly from numpy — no ID conversion
        cpd_counts_contracted = result['x_contracted'].sum(axis=1)
        rxn_counts_contracted = result['y_contracted'].sum(axis=1)
        cpd_counts_reexp = result['x_reexpanded'].sum(axis=1)
        rxn_counts_reexp = result['y_reexpanded'].sum(axis=1)

        # Verify counts match the length of ID lists from separate calls
        cpds_c, rxns_c = toy._run_contractions_rust(rxns_scope, cpds_scope, extinct_sets, seeds)
        for i in range(len(extinct_sets)):
            self.assertEqual(cpd_counts_contracted[i], len(cpds_c[i]))
            self.assertEqual(rxn_counts_contracted[i], len(rxns_c[i]))

    def test_seeds_preserved_in_fused_call(self):
        """Seeds survive through the fused contract_and_reexpand."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C"], ["D"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        # Kill rxn 0 (A + B -> C)
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 0]

        result = toy.contract_and_reexpand(seeds, rxns_scope, cpds_scope, [rxns_to_kill])

        # Seeds must be present in contracted output
        x_contracted = result['x_contracted'][0]
        seed_indices = [toy.cid_to_idx[s] for s in seeds if s in toy.cid_to_idx]
        for idx in seed_indices:
            self.assertEqual(x_contracted[idx], 1,
                             f"Seed at index {idx} missing from contracted output")


class TestContractAndReexpandKEGG(unittest.TestCase):
    """Tests on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
        cls.cpds_scope, cls.rxns_scope = cls.kegg.expand(cls.seeds)
        cls.kegg._ensure_dicts()
        cls.all_rxns = list(cls.kegg.rid_to_idx.keys())

    def test_kegg_matches_separate_calls(self):
        """Fused call matches separate contraction + expansion on KEGG."""
        random_seed(42)
        extinct_sets = [sample(self.all_rxns, len(self.all_rxns) // 10) for _ in range(5)]

        result = self.kegg.contract_and_reexpand(
            self.seeds, self.rxns_scope, self.cpds_scope, extinct_sets)

        # Compare contraction output
        cpds_c, rxns_c = self.kegg._run_contractions_rust(
            self.rxns_scope, self.cpds_scope, extinct_sets, self.seeds)
        for i in range(5):
            expected_x = self.kegg.initialize_metabolite_vector(cpds_c[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['x_contracted'][i], expected_x)

        # Compare re-expansion output
        cpds_r, rxns_r = self.kegg.run_expansions_batch(cpds_c, extinct_sets)
        for i in range(5):
            expected_x = self.kegg.initialize_metabolite_vector(cpds_r[i]).astype(np.uint8)
            np.testing.assert_array_equal(result['x_reexpanded'][i], expected_x)

    def test_kegg_all_seeds_preserved(self):
        """All seeds survive through fused contract_and_reexpand on KEGG."""
        random_seed(99)
        extinct_sets = [sample(self.all_rxns, len(self.all_rxns) // 5) for _ in range(3)]

        result = self.kegg.contract_and_reexpand(
            self.seeds, self.rxns_scope, self.cpds_scope, extinct_sets)

        seeds_in_scope = set(self.seeds) & set(self.cpds_scope)
        seed_indices = [self.kegg.cid_to_idx[s] for s in seeds_in_scope]

        for i in range(3):
            for idx in seed_indices:
                self.assertEqual(result['x_contracted'][i, idx], 1)

    def test_kegg_large_batch(self):
        """Large batch (100 extinction sets) runs without error."""
        random_seed(123)
        base_rxns = list(set(r[0] for r in self.all_rxns))
        extinct_sets = [sample(base_rxns, len(base_rxns) // 10) for _ in range(100)]

        result = self.kegg.contract_and_reexpand(
            self.seeds, self.rxns_scope, self.cpds_scope, extinct_sets)

        self.assertEqual(result['x_contracted'].shape[0], 100)
        self.assertEqual(result['x_reexpanded'].shape[0], 100)
        # Re-expansion should always produce >= contracted compounds
        for i in range(100):
            self.assertGreaterEqual(
                result['x_reexpanded'][i].sum(),
                result['x_contracted'][i].sum())


if __name__ == "__main__":
    unittest.main()
