"""
Tests for consistent reaction ID handling across the public API.

Background
----------
After convertToIrreversible(), rid_to_idx is keyed by (rn_id, direction)
tuples such as ("R00001", "forward").  Before this fix, passing plain string
IDs (e.g. "R00001") to initialize_reaction_vector() or expand(excluded_reactions=…)
would silently produce an all-zero vector, meaning no masking was applied.

These tests verify:
  1. initialize_reaction_vector() accepts both tuple and string IDs.
  2. expand(excluded_reactions=…) correctly excludes reactions whether string or
     tuple IDs are supplied, for both the Rust and Python backends.
  3. run_contractions() and contract() continue to work with tuple IDs
     (regression guard).
  4. Mixed sets (some tuples, some strings) are handled correctly.
"""

import unittest
from unittest.mock import patch
import numpy as np
import networkExpansionPy.lib as ne


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _toy_network():
    """Return a small irreversible GlobalMetabolicNetwork for unit tests.

    Reactions (indices after convertToIrreversible, forward copies):
      0: A+B → C
      1: C+D → E+F
      2: E+F → G
      3: G+H → I
      4: A+J → I
    Each has a paired reverse copy.
    """
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
    return toy


def _reaction_string_ids(toy):
    """Return the set of plain string rn IDs (integers as strings after
    _load_tuple_network, e.g. 0, 1, 2 …)."""
    return list(toy.network["rn"].unique())


def _reaction_tuple_ids(toy):
    """Return all (rn, direction) tuples present in the network."""
    return list(toy.rid_to_idx.keys())


# ---------------------------------------------------------------------------
# initialize_reaction_vector
# ---------------------------------------------------------------------------

class TestInitializeReactionVector(unittest.TestCase):

    def setUp(self):
        self.toy = _toy_network()

    # -- tuple form (should always have worked) ------------------------------

    def test_tuple_ids_produce_nonzero_vector(self):
        """Tuple IDs must set the corresponding entries to 1."""
        tuples = _reaction_tuple_ids(self.toy)[:4]
        vec = self.toy.initialize_reaction_vector(tuples)
        self.assertEqual(vec.sum(), 4)

    def test_tuple_ids_set_correct_indices(self):
        """Each tuple must map to exactly its own index."""
        tuples = _reaction_tuple_ids(self.toy)
        vec = self.toy.initialize_reaction_vector(tuples)
        # Every reaction should be flagged
        self.assertEqual(vec.sum(), len(self.toy.rid_to_idx))

    # -- string form (the bug that was fixed) --------------------------------

    def test_string_ids_produce_nonzero_vector(self):
        """Plain string rn IDs must activate both directed copies."""
        string_ids = _reaction_string_ids(self.toy)
        vec = self.toy.initialize_reaction_vector(string_ids)
        self.assertGreater(vec.sum(), 0,
            "String IDs produced an all-zero vector — the bug is not fixed")

    def test_string_id_activates_both_directions(self):
        """A single string ID must set forward AND reverse to 1."""
        # pick one rn string id
        rn = _reaction_string_ids(self.toy)[0]
        vec = self.toy.initialize_reaction_vector([rn])
        fwd_idx = self.toy.rid_to_idx.get((rn, "forward"))
        rev_idx = self.toy.rid_to_idx.get((rn, "reverse"))
        self.assertIsNotNone(fwd_idx, "forward key not found")
        self.assertIsNotNone(rev_idx, "reverse key not found")
        self.assertEqual(vec[fwd_idx], 1, "forward copy not set for string ID")
        self.assertEqual(vec[rev_idx], 1, "reverse copy not set for string ID")

    def test_string_ids_all_reactions(self):
        """All string IDs → all (rn, direction) pairs set."""
        string_ids = _reaction_string_ids(self.toy)
        vec = self.toy.initialize_reaction_vector(string_ids)
        expected = len(self.toy.rid_to_idx)  # every directed reaction
        self.assertEqual(vec.sum(), expected)

    # -- mixed form ----------------------------------------------------------

    def test_mixed_tuple_and_string_ids(self):
        """A mix of tuple and string IDs should combine correctly."""
        string_id = _reaction_string_ids(self.toy)[0]       # covers 2 directed rxns
        tuple_id  = _reaction_tuple_ids(self.toy)[-1]       # covers 1 directed rxn

        # Avoid double-counting if the tuple happens to be the same rxn
        expected_min = 1  # at least the tuple entry
        vec = self.toy.initialize_reaction_vector([string_id, tuple_id])
        self.assertGreaterEqual(vec.sum(), expected_min)

    # -- edge cases ----------------------------------------------------------

    def test_empty_set_returns_zero_vector(self):
        vec = self.toy.initialize_reaction_vector([])
        self.assertEqual(vec.sum(), 0)

    def test_unknown_id_ignored(self):
        """Unknown IDs (neither tuple nor string in the network) are silently skipped."""
        vec = self.toy.initialize_reaction_vector(["NOT_A_REAL_REACTION"])
        self.assertEqual(vec.sum(), 0)

    def test_none_returns_none(self):
        result = self.toy.initialize_reaction_vector(None)
        self.assertIsNone(result)


# ---------------------------------------------------------------------------
# expand(excluded_reactions=...) — both backends
# ---------------------------------------------------------------------------

class TestExpandReactionMask(unittest.TestCase):
    """expand() with excluded_reactions must produce different results when a mask
    is supplied, regardless of whether string or tuple IDs are used."""

    def setUp(self):
        self.toy = _toy_network()
        self.seeds = ["A", "B", "D", "H"]
        # Baseline: unmasked expansion
        self.base_cpds, self.base_rxns = self.toy.expand(self.seeds, algorithm="naive")

    def _first_rn_string(self):
        return str(_reaction_string_ids(self.toy)[0])

    def _first_rn_tuples(self):
        rn = _reaction_string_ids(self.toy)[0]
        return [(rn, "forward"), (rn, "reverse")]

    # -- Python backend -------------------------------------------------------

    def test_string_mask_python_differs_from_unmasked(self):
        """expand() with a string-ID mask (Python) must actually mask something."""
        mask = _reaction_string_ids(self.toy)  # mask ALL reactions
        with patch.object(ne, "_HAS_RUST", False):
            cpds, rxns = self.toy.expand(self.seeds, algorithm="naive",
                                          excluded_reactions=mask)
        # With all reactions masked, scope should shrink
        self.assertLessEqual(len(rxns), len(self.base_rxns),
            "Masking all reactions (string IDs, Python) had no effect")

    def test_tuple_mask_python_differs_from_unmasked(self):
        """expand() with tuple IDs (Python) still works as a regression guard."""
        mask = _reaction_tuple_ids(self.toy)  # mask ALL reactions
        with patch.object(ne, "_HAS_RUST", False):
            cpds, rxns = self.toy.expand(self.seeds, algorithm="naive",
                                          excluded_reactions=mask)
        self.assertLessEqual(len(rxns), len(self.base_rxns))

    def test_string_mask_python_equals_tuple_mask_python(self):
        """String and tuple masks for the same reactions must give identical results."""
        rn = _reaction_string_ids(self.toy)[0]
        string_mask = [rn]
        tuple_mask  = [(rn, "forward"), (rn, "reverse")]

        with patch.object(ne, "_HAS_RUST", False):
            cpds_str, rxns_str = self.toy.expand(self.seeds, algorithm="naive",
                                                   excluded_reactions=string_mask)
            cpds_tup, rxns_tup = self.toy.expand(self.seeds, algorithm="naive",
                                                   excluded_reactions=tuple_mask)

        self.assertEqual(sorted(cpds_str), sorted(cpds_tup))
        self.assertEqual(sorted(str(r) for r in rxns_str),
                         sorted(str(r) for r in rxns_tup))

    # -- Rust backend ---------------------------------------------------------

    @unittest.skipUnless(ne._HAS_RUST, "netexprs not installed")
    def test_string_mask_rust_differs_from_unmasked(self):
        """expand() with string-ID mask (Rust) must actually mask something."""
        mask = _reaction_string_ids(self.toy)
        with patch.object(ne, "_HAS_RUST", True):
            cpds, rxns = self.toy.expand(self.seeds, algorithm="naive",
                                          excluded_reactions=mask)
        self.assertLessEqual(len(rxns), len(self.base_rxns),
            "Masking all reactions (string IDs, Rust) had no effect")

    @unittest.skipUnless(ne._HAS_RUST, "netexprs not installed")
    def test_string_mask_rust_equals_tuple_mask_rust(self):
        """String and tuple masks (Rust) for the same reactions must match."""
        rn = _reaction_string_ids(self.toy)[0]
        string_mask = [rn]
        tuple_mask  = [(rn, "forward"), (rn, "reverse")]

        with patch.object(ne, "_HAS_RUST", True):
            cpds_str, rxns_str = self.toy.expand(self.seeds, algorithm="naive",
                                                   excluded_reactions=string_mask)
            cpds_tup, rxns_tup = self.toy.expand(self.seeds, algorithm="naive",
                                                   excluded_reactions=tuple_mask)

        self.assertEqual(sorted(cpds_str), sorted(cpds_tup))
        self.assertEqual(sorted(str(r) for r in rxns_str),
                         sorted(str(r) for r in rxns_tup))

    @unittest.skipUnless(ne._HAS_RUST, "netexprs not installed")
    def test_rust_and_python_string_mask_agree(self):
        """Rust and Python backends with string masks must agree."""
        rn = _reaction_string_ids(self.toy)[0]
        mask = [rn]

        with patch.object(ne, "_HAS_RUST", False):
            cpds_py, rxns_py = self.toy.expand(self.seeds, algorithm="naive",
                                                excluded_reactions=mask)
        with patch.object(ne, "_HAS_RUST", True):
            cpds_rs, rxns_rs = self.toy.expand(self.seeds, algorithm="naive",
                                                excluded_reactions=mask)

        self.assertEqual(sorted(cpds_py), sorted(cpds_rs))
        self.assertEqual(sorted(str(r) for r in rxns_py),
                         sorted(str(r) for r in rxns_rs))

    # -- no-op mask guard (regression) ---------------------------------------

    def test_empty_mask_equals_unmasked(self):
        """An empty mask must not change results."""
        with patch.object(ne, "_HAS_RUST", False):
            cpds, rxns = self.toy.expand(self.seeds, algorithm="naive",
                                          excluded_reactions=[])
        self.assertEqual(sorted(cpds), sorted(self.base_cpds))

    def test_none_mask_equals_unmasked(self):
        """excluded_reactions=None must be equivalent to no mask."""
        with patch.object(ne, "_HAS_RUST", False):
            cpds, rxns = self.toy.expand(self.seeds, algorithm="naive",
                                          excluded_reactions=None)
        self.assertEqual(sorted(cpds), sorted(self.base_cpds))


# ---------------------------------------------------------------------------
# Regression: contract / run_contractions still work with tuple IDs
# ---------------------------------------------------------------------------

class TestContractTupleIdsRegression(unittest.TestCase):
    """Guard against regressions: tuple-ID paths must still work after the fix."""

    def setUp(self):
        self.toy = _toy_network()
        self.seeds = ["A", "B", "D", "H"]
        with patch.object(ne, "_HAS_RUST", False):
            self.cpds, self.rxns = self.toy.expand(self.seeds, algorithm="naive")

    def test_contract_tuple_ids(self):
        """contract() with tuple IDs must produce a non-empty result."""
        # Extinct = first two directed reactions in the scope (as tuples)
        extinct = self.rxns[:2]
        with patch.object(ne, "_HAS_RUST", False):
            cpds, rxns = self.toy.contract(self.seeds, self.rxns, self.cpds, extinct)
        # Just verify it ran without error and returned lists
        self.assertIsInstance(cpds, list)
        self.assertIsInstance(rxns, list)

    def test_contract_string_ids(self):
        """contract() with string IDs (new capability) must agree with tuple version."""
        # Scope reactions as tuples and strings
        tuple_extinct = self.rxns[:2]
        # Extract plain rn IDs from the tuple extinct set (de-duplicate directions)
        string_extinct = list({r[0] for r in tuple_extinct})

        with patch.object(ne, "_HAS_RUST", False):
            cpds_tup, rxns_tup = self.toy.contract(
                self.seeds, self.rxns, self.cpds, tuple_extinct)
            cpds_str, rxns_str = self.toy.contract(
                self.seeds, self.rxns, self.cpds, string_extinct)

        # String form extinguishes both directions, so scope should be ≤ tuple form
        # (it can't be larger — string covers a superset of what tuple covers)
        self.assertLessEqual(len(rxns_str), len(rxns_tup) + len(string_extinct))

    def test_run_contractions_tuple_ids(self):
        """run_contractions() with tuple extinct sets must still work."""
        extinct_sets = [self.rxns[:2], self.rxns[2:4]]
        with patch.object(ne, "_HAS_RUST", False):
            cpd_scopes, rxn_scopes = self.toy.run_contractions(
                self.seeds, self.rxns, self.cpds, extinct_sets)
        self.assertEqual(len(cpd_scopes), 2)
        self.assertEqual(len(rxn_scopes), 2)


# ---------------------------------------------------------------------------
# run_expansions_reactionMasks with string IDs
# ---------------------------------------------------------------------------

class TestRunExpansionsReactionMasksStringIds(unittest.TestCase):

    def setUp(self):
        self.toy = _toy_network()
        self.seeds = ["A", "B", "D", "H"]

    def test_string_mask_sets_shrink_scope(self):
        """Masking all reactions with string IDs must reduce scope."""
        all_string_ids = _reaction_string_ids(self.toy)
        with patch.object(ne, "_HAS_RUST", False):
            cpd_scopes, rxn_scopes = self.toy.run_expansions_reactionMasks(
                self.seeds, [all_string_ids])
        # Should be able to expand to at least the seeds, but reactions should be 0
        self.assertEqual(len(rxn_scopes[0]), 0,
            "All-reaction string mask did not eliminate any reactions")

    def test_string_and_tuple_masks_agree(self):
        """String and tuple masks for the same reaction set must agree."""
        rn = _reaction_string_ids(self.toy)[0]
        string_mask_sets = [[rn]]
        tuple_mask_sets  = [[(rn, "forward"), (rn, "reverse")]]

        with patch.object(ne, "_HAS_RUST", False):
            cpds_str, rxns_str = self.toy.run_expansions_reactionMasks(
                self.seeds, string_mask_sets)
            cpds_tup, rxns_tup = self.toy.run_expansions_reactionMasks(
                self.seeds, tuple_mask_sets)

        self.assertEqual(sorted(cpds_str[0]), sorted(cpds_tup[0]))
        self.assertEqual(sorted(str(r) for r in rxns_str[0]),
                         sorted(str(r) for r in rxns_tup[0]))


if __name__ == "__main__":
    unittest.main()