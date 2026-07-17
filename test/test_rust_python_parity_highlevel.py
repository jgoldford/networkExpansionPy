"""
High-level parity tests for Rust vs Python network expansion backends.

These tests call the GlobalMetabolicNetwork class methods directly (e.g., expand(),
run_expansions(), run_contractions()) rather than the low-level array functions.
This ensures the full pipeline (ID conversion, array preparation, result conversion)
produces identical results between Rust and Python backends.

The low-level tests in test_rust_python_parity_lowlevel.py isolate the Rust
implementation from the Python wrapper logic. These high-level tests ensure
the integration layer doesn't introduce bugs (e.g., mask semantics being
inverted, ID lookups being wrong, etc.).
"""

import unittest
import numpy as np
from random import sample, seed as random_seed
import networkExpansionPy.lib as ne


class TestExpandHighLevel(unittest.TestCase):
    """Test single expansion through the high-level API."""

    def test_toy_network_expand_parity(self):
        """Compare _expand_rust vs _expand_python on toy network."""
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

        seedSet = ["A", "B", "D", "H"]

        cpds_rust, rxns_rust = toy._expand_rust(seedSet)
        cpds_py, rxns_py = toy._expand_python(seedSet, algorithm='naive')

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_toy_network_expand_minimal_seed(self):
        """Single compound seed that can't fire any reactions."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        seedSet = ["A"]  # Can't fire A+B->C without B

        cpds_rust, rxns_rust = toy._expand_rust(seedSet)
        cpds_py, rxns_py = toy._expand_python(seedSet, algorithm='naive')

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))
        # Should only contain the seed
        self.assertEqual(set(cpds_rust), {"A"})
        self.assertEqual(set(rxns_rust), set())

    def test_toy_network_expand_with_mask_exclude(self):
        """Test masked expansion where mask specifies reactions to EXCLUDE."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),       # rxn 0
            (["C", "D"], ["E", "F"]),  # rxn 1
            (["E", "F"], ["G"]),       # rxn 2
            (["G", "H"], ["I"]),       # rxn 3
            (["A", "J"], ["I"]),       # rxn 4
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seedSet = ["A", "B", "D", "H"]

        # Get all reaction IDs for rxn 1 (both forward and reverse)
        # This should block C+D -> E,F, limiting the expansion
        rxns_to_exclude = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 1]

        cpds_rust, rxns_rust = toy._expand_rust(seedSet, excluded_reactions=rxns_to_exclude)
        cpds_py, rxns_py = toy._expand_python(seedSet, algorithm='naive', excluded_reactions=rxns_to_exclude)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_toy_network_expand_empty_mask(self):
        """Empty mask should not exclude any reactions."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        seedSet = ["A", "B", "D"]

        # Empty list means exclude nothing
        cpds_rust_masked, rxns_rust_masked = toy._expand_rust(seedSet, excluded_reactions=[])
        cpds_rust_unmasked, rxns_rust_unmasked = toy._expand_rust(seedSet, excluded_reactions=None)
        cpds_py, rxns_py = toy._expand_python(seedSet, algorithm='naive', excluded_reactions=[])

        # All three should produce the same results
        self.assertEqual(set(cpds_rust_masked), set(cpds_py))
        self.assertEqual(set(rxns_rust_masked), set(rxns_py))
        self.assertEqual(set(cpds_rust_masked), set(cpds_rust_unmasked))
        self.assertEqual(set(rxns_rust_masked), set(rxns_rust_unmasked))


class TestExpandHighLevelLargeNetwork(unittest.TestCase):
    """Test single expansion on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_expand_parity(self):
        """Compare _expand_rust vs _expand_python on KEGG network."""
        cpds_rust, rxns_rust = self.kegg._expand_rust(self.seedSet)
        cpds_py, rxns_py = self.kegg._expand_python(self.seedSet, algorithm='naive')

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_kegg_expand_minimal_seed(self):
        """Single compound seed on KEGG."""
        seedSet = ["C00001"]  # Water

        cpds_rust, rxns_rust = self.kegg._expand_rust(seedSet)
        cpds_py, rxns_py = self.kegg._expand_python(seedSet, algorithm='naive')

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_kegg_expand_with_mask(self):
        """Masked expansion on KEGG network."""
        self.kegg._ensure_dicts()

        # Randomly select 10% of reactions to exclude
        all_rxns = list(self.kegg.rid_to_idx.keys())
        np.random.seed(42)
        n_exclude = len(all_rxns) // 10
        rxns_to_exclude = sample(all_rxns, n_exclude)

        cpds_rust, rxns_rust = self.kegg._expand_rust(self.seedSet, excluded_reactions=rxns_to_exclude)
        cpds_py, rxns_py = self.kegg._expand_python(self.seedSet, algorithm='naive', excluded_reactions=rxns_to_exclude)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))


class TestContractHighLevel(unittest.TestCase):
    """Test single contraction through the high-level API."""

    def test_toy_network_contract_parity(self):
        """Compare _contract_rust vs _contract_python on toy network."""
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

        seedSet = ["A", "B", "D", "H"]

        # First expand to get full scope
        cpds_scope, rxns_scope = toy._expand_python(seedSet, algorithm='naive')

        # Get some reactions to mark as extinct
        toy._ensure_dicts()
        all_rxns = list(toy.rid_to_idx.keys())
        extinct_rxns = all_rxns[:2]  # First two reactions

        cpds_rust, rxns_rust = toy._contract_rust(rxns_scope, cpds_scope, extinct_rxns)
        cpds_py, rxns_py = toy._contract_python(rxns_scope, cpds_scope, extinct_rxns)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_toy_network_contract_no_extinction(self):
        """Contract with no extinct reactions (should be no-op)."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C", "D"], ["E", "F"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()

        seedSet = ["A", "B", "D"]
        cpds_scope, rxns_scope = toy._expand_python(seedSet, algorithm='naive')

        # No extinctions
        cpds_rust, rxns_rust = toy._contract_rust(rxns_scope, cpds_scope, [])
        cpds_py, rxns_py = toy._contract_python(rxns_scope, cpds_scope, [])

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))


class TestContractHighLevelLargeNetwork(unittest.TestCase):
    """Test single contraction on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
        # Pre-compute expansion scope
        cls.cpds_scope, cls.rxns_scope = cls.kegg._expand_python(cls.seedSet, algorithm='naive')

    def test_kegg_contract_parity(self):
        """Compare _contract_rust vs _contract_python on KEGG."""
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())

        # Randomly select 5% of reactions to extinct
        np.random.seed(42)
        n_extinct = len(all_rxns) // 20
        extinct_rxns = sample(all_rxns, n_extinct)

        cpds_rust, rxns_rust = self.kegg._contract_rust(self.rxns_scope, self.cpds_scope, extinct_rxns)
        cpds_py, rxns_py = self.kegg._contract_python(self.rxns_scope, self.cpds_scope, extinct_rxns)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_kegg_contract_heavy_extinction(self):
        """Contract with many extinct reactions."""
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())

        # Randomly select 30% of reactions to extinct
        np.random.seed(123)
        n_extinct = int(len(all_rxns) * 0.3)
        extinct_rxns = sample(all_rxns, n_extinct)

        cpds_rust, rxns_rust = self.kegg._contract_rust(self.rxns_scope, self.cpds_scope, extinct_rxns)
        cpds_py, rxns_py = self.kegg._contract_python(self.rxns_scope, self.cpds_scope, extinct_rxns)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))


class TestRunExpansionsHighLevel(unittest.TestCase):
    """Test batch expansions through run_expansions()."""

    def test_toy_batch_expansions_parity(self):
        """Compare Rust vs Python batch expansions on toy network."""
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

        seedSets = [
            ["A", "B", "D", "H"],
            ["A", "B"],
            ["C", "D"],
            ["I"],  # reverse direction
            ["A", "B", "D", "H", "J"],  # full expansion
        ]

        # run_expansions uses Rust when available, compare against pure Python
        cpds_rust_list, rxns_rust_list = toy.run_expansions(seedSets, algorithm='naive')
        cpds_py_list, rxns_py_list = toy._run_expansions_python(seedSets, algorithm='naive')

        self.assertEqual(len(cpds_rust_list), len(cpds_py_list))
        for i in range(len(seedSets)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for seed set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for seed set {i}")


class TestRunExpansionsHighLevelLargeNetwork(unittest.TestCase):
    """Test batch expansions on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_kegg_batch_expansions_parity(self):
        """Compare Rust vs Python batch expansions on KEGG."""
        # Create varying seed sets
        random_seed(42)
        seedSets = [
            self.seedSet,
            self.seedSet[:10],
            self.seedSet[:5],
            ["C00001"],  # Water only
            self.seedSet + ["C00002", "C00003"],
        ]

        # run_expansions uses Rust when available, compare against pure Python
        cpds_rust_list, rxns_rust_list = self.kegg.run_expansions(seedSets, algorithm='naive')
        cpds_py_list, rxns_py_list = self.kegg._run_expansions_python(seedSets, algorithm='naive')

        self.assertEqual(len(cpds_rust_list), len(cpds_py_list))
        for i in range(len(seedSets)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for seed set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for seed set {i}")

    def test_kegg_batch_expansions_many_seeds(self):
        """Batch expansion with many different seed sets."""
        random_seed(123)
        n_batches = 20
        seedSets = []
        for _ in range(n_batches):
            n_seeds = np.random.randint(1, len(self.seedSet))
            seedSets.append(sample(self.seedSet, n_seeds))

        # run_expansions uses Rust when available, compare against pure Python
        cpds_rust_list, rxns_rust_list = self.kegg.run_expansions(seedSets, algorithm='naive')
        cpds_py_list, rxns_py_list = self.kegg._run_expansions_python(seedSets, algorithm='naive')

        for i in range(n_batches):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for seed set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for seed set {i}")


class TestRunContractionsHighLevel(unittest.TestCase):
    """Test batch contractions through run_contractions()."""

    def test_toy_batch_contractions_parity(self):
        """Compare Rust vs Python batch contractions on toy network."""
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

        seedSet = ["A", "B", "D", "H"]
        cpds_scope, rxns_scope = toy._expand_python(seedSet, algorithm='naive')

        toy._ensure_dicts()
        all_rxns = list(toy.rid_to_idx.keys())

        # Create varying extinction sets
        extinctSets = [
            [],
            all_rxns[:1],
            all_rxns[:3],
            all_rxns[::2],  # Every other reaction
            all_rxns,  # All reactions extinct
        ]

        cpds_rust_list, rxns_rust_list = toy._run_contractions_rust(rxns_scope, cpds_scope, extinctSets)
        cpds_py_list, rxns_py_list = toy._run_contractions_python(rxns_scope, cpds_scope, extinctSets)

        self.assertEqual(len(cpds_rust_list), len(cpds_py_list))
        for i in range(len(extinctSets)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for extinction set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for extinction set {i}")


class TestRunContractionsHighLevelLargeNetwork(unittest.TestCase):
    """Test batch contractions on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
        cls.cpds_scope, cls.rxns_scope = cls.kegg._expand_python(cls.seedSet, algorithm='naive')
        cls.kegg._ensure_dicts()
        cls.all_rxns = list(cls.kegg.rid_to_idx.keys())

    def test_kegg_batch_contractions_parity(self):
        """Compare Rust vs Python batch contractions on KEGG."""
        random_seed(42)
        n_batches = 10
        extinctSets = []
        for _ in range(n_batches):
            n_extinct = np.random.randint(0, len(self.all_rxns) // 10)
            extinctSets.append(sample(self.all_rxns, n_extinct))

        cpds_rust_list, rxns_rust_list = self.kegg._run_contractions_rust(
            self.rxns_scope, self.cpds_scope, extinctSets)
        cpds_py_list, rxns_py_list = self.kegg._run_contractions_python(
            self.rxns_scope, self.cpds_scope, extinctSets)

        for i in range(n_batches):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for extinction set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for extinction set {i}")

    def test_kegg_batch_contractions_varying_rates(self):
        """Batch contraction with varying extinction rates."""
        random_seed(456)
        extinction_rates = [0.01, 0.05, 0.1, 0.2, 0.3]
        extinctSets = []
        for rate in extinction_rates:
            n_extinct = int(len(self.all_rxns) * rate)
            extinctSets.append(sample(self.all_rxns, n_extinct))

        cpds_rust_list, rxns_rust_list = self.kegg._run_contractions_rust(
            self.rxns_scope, self.cpds_scope, extinctSets)
        cpds_py_list, rxns_py_list = self.kegg._run_contractions_python(
            self.rxns_scope, self.cpds_scope, extinctSets)

        for i in range(len(extinction_rates)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for extinction rate {extinction_rates[i]}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for extinction rate {extinction_rates[i]}")


class TestRunExpansionsReactionMasksHighLevel(unittest.TestCase):
    """Test batch masked expansions through run_expansions_reactionMasks()."""

    def test_toy_batch_masked_expansions_parity(self):
        """Compare Rust vs Python batch masked expansions on toy network."""
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

        seedSet = ["A", "B", "D", "H"]
        all_rxns = list(toy.rid_to_idx.keys())

        # Create varying mask sets (reactions to REMOVE)
        maskedSets = [
            [],  # No masking
            all_rxns[:2],  # First two reactions
            all_rxns[::2],  # Every other reaction
            [r for r in all_rxns if r[0] == 0],  # All variants of reaction 0
        ]

        cpds_rust_list, rxns_rust_list = toy.run_expansions_batch(seedSet, maskedSets)
        cpds_py_list, rxns_py_list = toy._run_expansions_reactionMasks_python(seedSet, maskedSets)

        self.assertEqual(len(cpds_rust_list), len(cpds_py_list))
        for i in range(len(maskedSets)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for mask set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for mask set {i}")


class TestRunExpansionsReactionMasksHighLevelLargeNetwork(unittest.TestCase):
    """Test batch masked expansions on real KEGG network."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
        cls.kegg._ensure_dicts()
        cls.all_rxns = list(cls.kegg.rid_to_idx.keys())
        # Get unique reaction base names (without direction)
        cls.network_rxns = list(set(r[0] for r in cls.all_rxns))

    def test_kegg_batch_masked_expansions_parity(self):
        """Compare Rust vs Python batch masked expansions on KEGG."""
        random_seed(42)
        n_batches = 10
        maskedSets = []
        for _ in range(n_batches):
            # Sample base reaction names, then get all their direction variants
            n_remove = np.random.randint(0, len(self.network_rxns) // 10)
            rxns_removed_base = sample(self.network_rxns, n_remove)
            rxns_removed = [r for r in self.all_rxns if r[0] in rxns_removed_base]
            maskedSets.append(rxns_removed)

        cpds_rust_list, rxns_rust_list = self.kegg.run_expansions_batch(
            self.seedSet, maskedSets)
        cpds_py_list, rxns_py_list = self.kegg._run_expansions_reactionMasks_python(
            self.seedSet, maskedSets)

        for i in range(n_batches):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for mask set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for mask set {i}")

    def test_kegg_batch_masked_expansions_varying_fractions(self):
        """Batch masked expansion with varying removal fractions (like user's code)."""
        random_seed(789)

        # Similar to user's sampling code
        fractions = [0.01, 0.05, 0.1, 0.25, 0.5]
        reps = 3

        maskedSets = []
        for frac in fractions:
            sample_size = int(frac * len(self.network_rxns))
            for _ in range(reps):
                rxns_removed_base = sample(self.network_rxns, sample_size)
                rxns_removed = [r for r in self.all_rxns if r[0] in rxns_removed_base]
                maskedSets.append(rxns_removed)

        cpds_rust_list, rxns_rust_list = self.kegg.run_expansions_batch(
            self.seedSet, maskedSets)
        cpds_py_list, rxns_py_list = self.kegg._run_expansions_reactionMasks_python(
            self.seedSet, maskedSets)

        for i in range(len(maskedSets)):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for mask set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for mask set {i}")

    def test_kegg_batch_masked_expansions_large_batch(self):
        """Large batch of masked expansions."""
        random_seed(999)
        n_batches = 10

        maskedSets = []
        for _ in range(n_batches):
            n_remove = np.random.randint(0, len(self.network_rxns) // 5)
            rxns_removed_base = sample(self.network_rxns, n_remove)
            rxns_removed = [r for r in self.all_rxns if r[0] in rxns_removed_base]
            maskedSets.append(rxns_removed)

        cpds_rust_list, rxns_rust_list = self.kegg.run_expansions_batch(
            self.seedSet, maskedSets)
        cpds_py_list, rxns_py_list = self.kegg._run_expansions_reactionMasks_python(
            self.seedSet, maskedSets)

        for i in range(n_batches):
            self.assertEqual(set(cpds_rust_list[i]), set(cpds_py_list[i]),
                             f"Compound mismatch for mask set {i}")
            self.assertEqual(set(rxns_rust_list[i]), set(rxns_py_list[i]),
                             f"Reaction mismatch for mask set {i}")


class TestParallelMethodsHighLevel(unittest.TestCase):
    """Test that parallel methods produce same results as non-parallel."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
        cls.kegg._ensure_dicts()
        cls.all_rxns = list(cls.kegg.rid_to_idx.keys())
        cls.network_rxns = list(set(r[0] for r in cls.all_rxns))

    def test_run_expansions_parallel_parity(self):
        """run_expansions_parallel should match run_expansions."""
        random_seed(42)
        seedSets = [
            self.seedSet,
            self.seedSet[:10],
            sample(self.seedSet, 20),
        ]

        cpds_normal, rxns_normal = self.kegg.run_expansions(seedSets)
        cpds_parallel, rxns_parallel = self.kegg.run_expansions_parallel(seedSets)

        for i in range(len(seedSets)):
            self.assertEqual(set(cpds_normal[i]), set(cpds_parallel[i]),
                             f"Compound mismatch for seed set {i}")
            self.assertEqual(set(rxns_normal[i]), set(rxns_parallel[i]),
                             f"Reaction mismatch for seed set {i}")

    def test_run_expansions_reactionMasks_parallel_parity(self):
        """run_expansions_reactionMasks_parallel should match run_expansions_reactionMasks."""
        random_seed(123)
        n_batches = 10
        maskedSets = []
        for _ in range(n_batches):
            n_remove = np.random.randint(0, len(self.network_rxns) // 10)
            rxns_removed_base = sample(self.network_rxns, n_remove)
            rxns_removed = [r for r in self.all_rxns if r[0] in rxns_removed_base]
            maskedSets.append(rxns_removed)

        cpds_normal, rxns_normal = self.kegg.run_expansions_reactionMasks(self.seedSet, maskedSets)
        cpds_parallel, rxns_parallel = self.kegg.run_expansions_reactionMasks_parallel(
            self.seedSet, maskedSets)

        for i in range(n_batches):
            self.assertEqual(set(cpds_normal[i]), set(cpds_parallel[i]),
                             f"Compound mismatch for mask set {i}")
            self.assertEqual(set(rxns_normal[i]), set(rxns_parallel[i]),
                             f"Reaction mismatch for mask set {i}")


class TestPublicAPIDispatch(unittest.TestCase):
    """Test that public API methods correctly dispatch to Rust when available."""

    @classmethod
    def setUpClass(cls):
        cls.kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        cls.seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

    def test_expand_dispatches_correctly(self):
        """expand() should use Rust for naive algorithm and produce correct results."""
        # expand() with naive should dispatch to Rust
        cpds_api, rxns_api = self.kegg.expand(self.seedSet, algorithm='naive')
        cpds_py, rxns_py = self.kegg._expand_python(self.seedSet, algorithm='naive')

        self.assertEqual(set(cpds_api), set(cpds_py))
        self.assertEqual(set(rxns_api), set(rxns_py))

    def test_contract_dispatches_correctly(self):
        """contract() should use Rust and produce correct results."""
        cpds_scope, rxns_scope = self.kegg._expand_python(self.seedSet, algorithm='naive')
        self.kegg._ensure_dicts()
        all_rxns = list(self.kegg.rid_to_idx.keys())
        extinct_rxns = sample(all_rxns, len(all_rxns) // 20)

        cpds_api, rxns_api = self.kegg.contract(self.seedSet, rxns_scope, cpds_scope, extinct_rxns)
        cpds_py, rxns_py = self.kegg._contract_python(rxns_scope, cpds_scope, extinct_rxns, self.seedSet)

        self.assertEqual(set(cpds_api), set(cpds_py))
        self.assertEqual(set(rxns_api), set(rxns_py))


if __name__ == "__main__":
    unittest.main()