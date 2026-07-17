"""
Tests for seed preservation during network contraction.

Seeds are compounds that enter the network as initial conditions — they aren't
produced by any reaction. Previously, contraction incorrectly dropped them because
the recomputation step `x_active = (P @ y_active) > 0` only retains compounds
that are *products* of surviving reactions.

These tests verify:
1. Ground truth: seeds survive contraction even when no reaction produces them
2. Ground truth: seed preservation prevents cascading loss of dependent reactions
3. Parity: Rust and Python backends produce identical results with seeds
4. Backward compat: passing no seeds matches the old (buggy) behavior
"""

import unittest
import numpy as np
from scipy.sparse import csr_matrix
import networkExpansionPy.lib as ne


class TestSeedPreservationGroundTruth(unittest.TestCase):
    """Ground truth tests — verifies correct behavior, not just parity."""

    def test_seeds_survive_when_no_reaction_produces_them(self):
        """Seeds must remain in compound scope even if no surviving reaction produces them.

        Network: A + B -> C, C -> D
        Seeds: [A, B]
        Expand: scope = {A, B, C, D}
        Extinct: kill reaction 0 (A + B -> C, both directions)

        Without seed preservation: A and B vanish (no reaction produces them),
        which cascades to lose everything.
        With seed preservation: A and B remain, but C and D are lost since
        the reaction producing them is extinct.
        """
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),  # rxn 0
            (["C"], ["D"]),       # rxn 1
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B"]

        # Expand to get full scope
        cpds_scope, rxns_scope = toy.expand(seeds)
        self.assertIn("A", cpds_scope)
        self.assertIn("B", cpds_scope)
        self.assertIn("C", cpds_scope)
        self.assertIn("D", cpds_scope)

        # Kill reaction 0 (A + B -> C) in both directions
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 0]
        self.assertTrue(len(rxns_to_kill) > 0)

        # Contract WITH seeds
        cpds_contracted, rxns_contracted = toy.contract(
            seeds, rxns_scope, cpds_scope, rxns_to_kill
        )

        # Seeds must survive
        self.assertIn("A", cpds_contracted)
        self.assertIn("B", cpds_contracted)

    def test_seeds_prevent_cascade_extinction(self):
        """Preserving seeds prevents cascading loss of reactions that need them.

        Network: A -> B, B -> C, C -> D
        Seeds: [A]
        Expand: scope = {A, B, C, D}
        Extinct: kill rxn 0 AND rxn 1 (all directions)

        Surviving: rxn 2 fwd (C->D), rxn 2 rev (D->C)
        Without seeds: A is not produced by any surviving reaction -> lost.
            B is only produced by rxn 1 rev (extinct) -> lost.
            C and D survive (produced by rxn 2 rev and rxn 2 fwd respectively).
        With seeds: A is preserved -> still just {A, C, D} since B has no
            surviving production path. But A is the key difference.
        """
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A"], ["B"]),  # rxn 0
            (["B"], ["C"]),  # rxn 1
            (["C"], ["D"]),  # rxn 2
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A"]
        cpds_scope, rxns_scope = toy.expand(seeds)
        self.assertEqual(set(cpds_scope), {"A", "B", "C", "D"})

        # Kill rxn 0 and rxn 1 (both directions each)
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] in (0, 1)]

        # With seeds
        cpds_with_seeds, _ = toy.contract(seeds, rxns_scope, cpds_scope, rxns_to_kill)

        # Without seeds
        cpds_without_seeds, _ = toy._contract_rust(
            rxns_scope, cpds_scope, rxns_to_kill, seedSet=None)

        # Seed A must survive with seed preservation
        self.assertIn("A", cpds_with_seeds)
        # Without seeds, A is lost
        self.assertNotIn("A", cpds_without_seeds)
        # C and D survive in both cases (produced by rxn 2 fwd/rev)
        self.assertIn("C", cpds_with_seeds)
        self.assertIn("D", cpds_with_seeds)
        self.assertIn("C", cpds_without_seeds)
        self.assertIn("D", cpds_without_seeds)

    def test_seed_preservation_enables_downstream_survival(self):
        """When seeds are preserved, reactions depending on them can still fire.

        Network: A -> C, B -> C, C -> D
        Seeds: [A, B]
        Expand: scope = {A, B, C, D}
        Extinct: kill rxn 0 (A -> C)

        With seed preservation: A and B remain. Rxn 1 (B -> C) still fires
        because B is a seed. C survives, so rxn 2 (C -> D) fires. D survives.
        """
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A"], ["C"]),  # rxn 0 — will be killed
            (["B"], ["C"]),  # rxn 1 — survives, needs B (seed)
            (["C"], ["D"]),  # rxn 2 — survives if C survives
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        # Kill rxn 0 (A -> C)
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 0]

        cpds_contracted, rxns_contracted = toy.contract(
            seeds, rxns_scope, cpds_scope, rxns_to_kill
        )

        # Seeds survive
        self.assertIn("A", cpds_contracted)
        self.assertIn("B", cpds_contracted)
        # B -> C still works, so C survives
        self.assertIn("C", cpds_contracted)
        # C -> D still works, so D survives
        self.assertIn("D", cpds_contracted)

    def test_without_seeds_old_behavior_drops_nonproducts(self):
        """Without seed preservation, compounds not produced by any reaction vanish.

        This verifies the OLD (buggy) behavior still occurs when seedSet=None
        is passed internally, for backward compatibility of low-level functions.
        """
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),  # rxn 0
            (["C"], ["D"]),       # rxn 1
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        # Kill rxn 0 (A + B -> C)
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 0]

        # Call _contract_rust directly with seedSet=None (old behavior)
        cpds_no_seeds, rxns_no_seeds = toy._contract_rust(
            rxns_scope, cpds_scope, rxns_to_kill, seedSet=None
        )

        # Without seed preservation, A and B should vanish
        self.assertNotIn("A", cpds_no_seeds)
        self.assertNotIn("B", cpds_no_seeds)

    def test_all_seeds_preserved_kegg(self):
        """On real KEGG network, all seed compounds survive any contraction."""
        kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

        cpds_scope, rxns_scope = kegg.expand(seeds)

        # Randomly extinct 10% of reactions
        kegg._ensure_dicts()
        all_rxns = list(kegg.rid_to_idx.keys())
        np.random.seed(42)
        n_extinct = len(all_rxns) // 10
        extinct_rxns = list(np.random.choice(len(all_rxns), n_extinct, replace=False))
        extinct_rxns = [all_rxns[i] for i in extinct_rxns]

        cpds_contracted, rxns_contracted = kegg.contract(
            seeds, rxns_scope, cpds_scope, extinct_rxns
        )

        # Every seed that was in scope must survive
        seeds_in_scope = set(seeds) & set(cpds_scope)
        for seed in seeds_in_scope:
            self.assertIn(seed, cpds_contracted,
                          f"Seed {seed} was lost during contraction")

    def test_batch_contractions_preserve_seeds_kegg(self):
        """Batch contractions on KEGG preserve all seeds across all replicates."""
        kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

        cpds_scope, rxns_scope = kegg.expand(seeds)

        kegg._ensure_dicts()
        all_rxns = list(kegg.rid_to_idx.keys())
        np.random.seed(123)

        # 10 different extinction sets at varying rates
        extinct_sets = []
        for rate in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]:
            n_extinct = int(len(all_rxns) * rate)
            indices = np.random.choice(len(all_rxns), n_extinct, replace=False)
            extinct_sets.append([all_rxns[i] for i in indices])

        cpds_list, rxns_list = kegg.run_contractions(
            seeds, rxns_scope, cpds_scope, extinct_sets
        )

        seeds_in_scope = set(seeds) & set(cpds_scope)
        for i, cpds in enumerate(cpds_list):
            for seed in seeds_in_scope:
                self.assertIn(seed, cpds,
                              f"Seed {seed} lost in batch {i} (rate={0.05 + i*0.05:.2f})")


class TestSeedPreservationParity(unittest.TestCase):
    """Verify Rust and Python backends produce identical results with seeds."""

    def test_parity_toy_network(self):
        """Rust and Python match on toy network with seed preservation."""
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

        # Kill some reactions
        all_rxns = list(toy.rid_to_idx.keys())
        extinct = all_rxns[:4]

        cpds_rust, rxns_rust = toy._contract_rust(rxns_scope, cpds_scope, extinct, seedSet=seeds)
        cpds_py, rxns_py = toy._contract_python(rxns_scope, cpds_scope, extinct, seedSet=seeds)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_parity_kegg_network(self):
        """Rust and Python match on real KEGG network with seed preservation."""
        kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

        cpds_scope, rxns_scope = kegg.expand(seeds)

        kegg._ensure_dicts()
        all_rxns = list(kegg.rid_to_idx.keys())
        np.random.seed(99)
        n_extinct = len(all_rxns) // 10
        indices = np.random.choice(len(all_rxns), n_extinct, replace=False)
        extinct = [all_rxns[i] for i in indices]

        cpds_rust, rxns_rust = kegg._contract_rust(rxns_scope, cpds_scope, extinct, seedSet=seeds)
        cpds_py, rxns_py = kegg._contract_python(rxns_scope, cpds_scope, extinct, seedSet=seeds)

        self.assertEqual(set(cpds_rust), set(cpds_py))
        self.assertEqual(set(rxns_rust), set(rxns_py))

    def test_parity_batch_kegg(self):
        """Batch contraction parity with seeds on KEGG."""
        kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
        seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

        cpds_scope, rxns_scope = kegg.expand(seeds)

        kegg._ensure_dicts()
        all_rxns = list(kegg.rid_to_idx.keys())
        np.random.seed(77)

        extinct_sets = []
        for _ in range(5):
            n_extinct = np.random.randint(100, len(all_rxns) // 5)
            indices = np.random.choice(len(all_rxns), n_extinct, replace=False)
            extinct_sets.append([all_rxns[i] for i in indices])

        cpds_rust, rxns_rust = kegg._run_contractions_rust(
            rxns_scope, cpds_scope, extinct_sets, seedSet=seeds)
        cpds_py, rxns_py = kegg._run_contractions_python(
            rxns_scope, cpds_scope, extinct_sets, seedSet=seeds)

        for i in range(len(extinct_sets)):
            self.assertEqual(set(cpds_rust[i]), set(cpds_py[i]),
                             f"Compound mismatch at batch {i}")
            self.assertEqual(set(rxns_rust[i]), set(rxns_py[i]),
                             f"Reaction mismatch at batch {i}")


class TestSeedPreservationChangedBehavior(unittest.TestCase):
    """Tests that explicitly demonstrate the behavior change.

    These tests would FAIL with the old implementation (no seed preservation)
    and PASS with the new implementation.
    """

    def test_contraction_retains_more_compounds_than_without_seeds(self):
        """With seeds, contraction must retain at least as many compounds as without."""
        toy = ne.GlobalMetabolicNetwork("dev")
        rxns = [
            (["A", "B"], ["C"]),
            (["C"], ["D"]),
            (["D", "A"], ["E"]),
        ]
        toy.network = ne._load_tuple_network(rxns)
        toy.convertToIrreversible()
        toy._ensure_dicts()

        seeds = ["A", "B"]
        cpds_scope, rxns_scope = toy.expand(seeds)

        # Kill rxn 0 (A + B -> C)
        rxns_to_kill = [rid for rid in toy.rid_to_idx.keys() if rid[0] == 0]

        # With seeds
        cpds_with_seeds, _ = toy._contract_rust(
            rxns_scope, cpds_scope, rxns_to_kill, seedSet=seeds)

        # Without seeds
        cpds_without_seeds, _ = toy._contract_rust(
            rxns_scope, cpds_scope, rxns_to_kill, seedSet=None)

        # With seeds must retain strictly more
        self.assertGreater(len(cpds_with_seeds), len(cpds_without_seeds))
        # Specifically, seeds must be the difference
        self.assertTrue(set(seeds).issubset(set(cpds_with_seeds)))
        self.assertFalse(set(seeds).issubset(set(cpds_without_seeds)))

    def test_total_extinction_still_preserves_seeds(self):
        """Even if ALL reactions are extinct, seeds survive."""
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

        # Kill ALL reactions (both forward and reverse directions)
        all_rxns = list(toy.rid_to_idx.keys())
        self.assertGreater(len(all_rxns), 0)

        # Contract with seeds — pass all reactions in scope as extinct
        cpds_contracted, rxns_contracted = toy.contract(
            seeds, rxns_scope, cpds_scope, rxns_scope
        )

        # No reactions survive
        self.assertEqual(len(rxns_contracted), 0)
        # But seeds do!
        self.assertIn("A", cpds_contracted)
        self.assertIn("B", cpds_contracted)
        # Non-seed compounds are gone (no reactions to produce them)
        self.assertNotIn("C", cpds_contracted)
        self.assertNotIn("D", cpds_contracted)


if __name__ == "__main__":
    unittest.main()
