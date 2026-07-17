"""
Test suite to verify Rust/Python dispatch logic in networkExpansionPy.
"""

import unittest
from unittest.mock import patch, MagicMock
import numpy as np

# We'll patch at the module level
import networkExpansionPy.lib as ne


class TestRustDispatch(unittest.TestCase):
    """Test that methods correctly dispatch to Rust vs Python backends."""
    
    @classmethod
    def setUpClass(cls):
        """Load network once for all tests."""
        cls.kegg = ne.GlobalMetabolicNetwork()
        cls.kegg.pruneInconsistentReactions()
        cls.kegg.convertToIrreversible()
        cls.kegg._ensure_rust_ready()
        
        # Get some test data
        seeds_df = ne.pd.read_csv(ne.asset_path + "/compounds/seeds.Goldford2022.csv")
        cls.seeds = list(set(seeds_df["ID"]))[:50]
        
        # Run one expansion to get scope for contraction tests
        cls.compound_scope, cls.reaction_scope = cls.kegg.expand(cls.seeds)
    
    def test_expand_dispatches_to_rust_when_available(self):
        """expand() should call _expand_rust when _HAS_RUST is True."""
        with patch.object(self.kegg, '_expand_rust', wraps=self.kegg._expand_rust) as mock_rust, \
             patch.object(self.kegg, '_expand_python', wraps=self.kegg._expand_python) as mock_python:
            
            # Force Rust available
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.expand(self.seeds, algorithm='naive')
            
            mock_rust.assert_called_once()
            mock_python.assert_not_called()
    
    def test_expand_dispatches_to_python_when_rust_unavailable(self):
        """expand() should call _expand_python when _HAS_RUST is False."""
        with patch.object(self.kegg, '_expand_rust', wraps=self.kegg._expand_rust) as mock_rust, \
             patch.object(self.kegg, '_expand_python', wraps=self.kegg._expand_python) as mock_python:
            
            # Force Rust unavailable
            with patch.object(ne, '_HAS_RUST', False):
                self.kegg.expand(self.seeds, algorithm='naive')
            
            mock_rust.assert_not_called()
            mock_python.assert_called_once()
    
    def test_expand_dispatches_to_rust_for_trace_algorithm(self):
        """expand() should use _expand_trace_rust for 'trace' algorithm with Rust available."""
        with patch.object(self.kegg, '_expand_trace_rust', wraps=self.kegg._expand_trace_rust) as mock_trace, \
             patch.object(self.kegg, '_expand_python', wraps=self.kegg._expand_python) as mock_python:

            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.expand(self.seeds, algorithm='trace')

            mock_trace.assert_called_once()
            mock_python.assert_not_called()
    
    def test_contract_dispatches_to_rust_when_available(self):
        """contract() should call _contract_rust when _HAS_RUST is True."""
        extinct = self.reaction_scope[:10]
        
        with patch.object(self.kegg, '_contract_rust', wraps=self.kegg._contract_rust) as mock_rust, \
             patch.object(self.kegg, '_contract_python', wraps=self.kegg._contract_python) as mock_python:
            
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.contract(self.seeds, self.reaction_scope, self.compound_scope, extinct)
            
            mock_rust.assert_called_once()
            mock_python.assert_not_called()
    
    def test_contract_dispatches_to_python_when_rust_unavailable(self):
        """contract() should call _contract_python when _HAS_RUST is False."""
        extinct = self.reaction_scope[:10]
        
        with patch.object(self.kegg, '_contract_rust', wraps=self.kegg._contract_rust) as mock_rust, \
             patch.object(self.kegg, '_contract_python', wraps=self.kegg._contract_python) as mock_python:
            
            with patch.object(ne, '_HAS_RUST', False):
                self.kegg.contract(self.seeds, self.reaction_scope, self.compound_scope, extinct)
            
            mock_rust.assert_not_called()
            mock_python.assert_called_once()
    
    def test_run_contractions_dispatches_to_rust_when_available(self):
        """run_contractions() should call _run_contractions_rust when _HAS_RUST is True."""
        extinct_sets = [self.reaction_scope[:10], self.reaction_scope[10:20]]
        
        with patch.object(self.kegg, '_run_contractions_rust', wraps=self.kegg._run_contractions_rust) as mock_rust, \
             patch.object(self.kegg, '_run_contractions_python', wraps=self.kegg._run_contractions_python) as mock_python:
            
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.run_contractions(self.seeds, self.reaction_scope, self.compound_scope, extinct_sets)
            
            mock_rust.assert_called_once()
            mock_python.assert_not_called()
    
    def test_run_expansions_reactionMasks_dispatches_to_rust_when_available(self):
        """run_expansions_reactionMasks() should call Rust backend when available."""
        mask_sets = [self.reaction_scope[:100], self.reaction_scope[:200]]
        
        with patch.object(self.kegg, 'run_expansions_batch', 
                          wraps=self.kegg.run_expansions_batch) as mock_rust, \
             patch.object(self.kegg, '_run_expansions_reactionMasks_python', 
                          wraps=self.kegg._run_expansions_reactionMasks_python) as mock_python:
            
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.run_expansions_reactionMasks(self.seeds, mask_sets)
            
            mock_rust.assert_called_once()
            mock_python.assert_not_called()
    
    def test_parallel_methods_delegate_to_rust_when_available(self):
        """_parallel methods should delegate to non-parallel (Rust handles threading)."""
        mask_sets = [self.reaction_scope[:100]]
        
        with patch.object(self.kegg, 'run_expansions_reactionMasks', 
                          wraps=self.kegg.run_expansions_reactionMasks) as mock_main:
            
            with patch.object(ne, '_HAS_RUST', True):
                self.kegg.run_expansions_reactionMasks_parallel(self.seeds, mask_sets)
            
            # Should delegate to the main method
            mock_main.assert_called_once()


class TestCacheInvalidation(unittest.TestCase):
    """Test that cache is properly invalidated when network changes."""
    
    def setUp(self):
        """Fresh network for each test."""
        self.kegg = ne.GlobalMetabolicNetwork()
        self.kegg.pruneInconsistentReactions()
        self.kegg.convertToIrreversible()
    
    def test_cache_built_on_first_use(self):
        """Cache should be None initially, then built on first _ensure_rust_ready."""
        self.assertIsNone(self.kegg._rust_arrays)
        self.kegg._ensure_rust_ready()
        self.assertIsNotNone(self.kegg._rust_arrays)
    
    def test_cache_invalidated_on_network_assignment(self):
        """Assigning to .network should invalidate cache."""
        self.kegg._ensure_rust_ready()
        self.assertIsNotNone(self.kegg._rust_arrays)
        
        # Direct assignment should trigger invalidation
        self.kegg.network = self.kegg.network.head(1000)
        
        self.assertIsNone(self.kegg._rust_arrays)
    
    def test_cache_invalidated_by_prune_methods(self):
        """Pruning methods should invalidate cache."""
        self.kegg._ensure_rust_ready()
        self.assertIsNotNone(self.kegg._rust_arrays)
        
        self.kegg.subnetwork(list(self.kegg.network.rn.unique())[:100])
        
        self.assertIsNone(self.kegg._rust_arrays)
    
    def test_cache_rebuilt_after_invalidation(self):
        """Cache should be rebuilt on next use after invalidation."""
        self.kegg._ensure_rust_ready()
        original_n_reactions = self.kegg._rust_arrays['n_reactions']
        
        # Invalidate by pruning
        self.kegg.subnetwork(list(self.kegg.network.rn.unique())[:100])
        self.assertIsNone(self.kegg._rust_arrays)
        
        # Rebuild
        self.kegg._ensure_rust_ready()
        self.assertIsNotNone(self.kegg._rust_arrays)
        
        # Should have different size now
        self.assertNotEqual(self.kegg._rust_arrays['n_reactions'], original_n_reactions)


if __name__ == '__main__':
    unittest.main()