# networkExpansionPy

Python package to construct biosphere-level metabolic networks and run network expansion algorithms. This package contains functions to prune biochemical reactions based on thermodynamic constraints and simulate the metabolic scope of different chemical environments.

## Overview

The package is built from algorithms described in the following papers:

**Goldford J.E. et al**, *Remnants of an ancient metabolism without phosphate.* Cell 168, 1–9, March 9, 2017

**Goldford J.E. et al**, *Environmental boundary conditions for the origin of life converge to an organo-sulfur metabolism.* Nature Ecol Evo 3,12 1715-1724, November 11, 2019

## Installation

### Basic Installation

Clone the repository and install using pip in a conda or virtual environment:

```bash
git clone https://github.com/jgoldford/networkExpansionPy.git
cd networkExpansionPy
pip install -e .
```

### Rust Acceleration (Recommended)

For **50-1000x speedup** on batch operations, install the optional Rust backend:

```bash
pip install netexprs
```

The package automatically detects and uses the Rust backend when available, with seamless fallback to pure Python if not installed. No code changes required!

#### Performance Comparison

Benchmarks on KEGG metabolic network (8,743 compounds, 18,940 reactions) using an 8-core system:

| Operation | Python | Rust | Speedup |
|-----------|--------|------|---------|
| Single Expansion | 281 ms | 43 ms | **6.5x** |
| Single Masked Expansion | 2.2 s | 10 ms | **214x** |
| Batch Masked Expansion (n=10) | 34.2 s | 43 ms | **802x** |
| Batch Multi-Seed Expansion (n=10) | 4.1 s | 46 ms | **89x** |
| Single Trace | 316 ms | 45 ms | **7.1x** |
| Batch Trace (n=10) | 4.5 s | 48 ms | **94x** |
| Single Contraction | 34 ms | 3.3 ms | **10x** |
| Batch Contraction (n=10) | 4.2 s | 4 ms | **1034x** |

The Rust backend provides massive speedups for batch operations through parallel processing, making large-scale robustness analyses and parameter sweeps practical.

## Quick Start

### Basic Network Expansion

```python
import networkExpansionPy.lib as ne

# Load metabolic network
kegg = ne.GlobalMetabolicNetwork()
kegg.pruneInconsistentReactions()
kegg.convertToIrreversible()

# Define seed compounds
seeds = ["C00001", "C00002", "C00003"]  # Example compound IDs

# Run expansion
compound_scope, reaction_scope = kegg.expand(seeds)

print(f"Starting from {len(seeds)} compounds")
print(f"Expanded to {len(compound_scope)} compounds")
print(f"Using {len(reaction_scope)} reactions")
```

### Batch Masked Expansions

Test network robustness by removing different reaction sets:

```python
# Define multiple reaction masks (reactions to remove)
masked_reaction_sets = [
    ["R00001", "R00002"],  # Remove these reactions
    ["R00003", "R00004"],  # Or these
    # ... many more ...
]

# Run all expansions in parallel (automatic with Rust backend)
compound_scopes, reaction_scopes = kegg.run_expansions_reactionMasks_parallel(
    seeds, 
    masked_reaction_sets
)

# Analyze results
for i, (cpds, rxns) in enumerate(zip(compound_scopes, reaction_scopes)):
    print(f"Mask {i}: {len(cpds)} compounds, {len(rxns)} reactions")
```

### Network Contraction

Remove reactions and contract the network:

```python
# Start with a full expansion
compound_scope, reaction_scope = kegg.expand(seeds)

# Define reactions to extinct
extinct_reactions = ["R00001", "R00002"]

# Contract the network
remaining_compounds, remaining_reactions = kegg.contract(
    seeds,
    reaction_scope,
    compound_scope, 
    extinct_reactions
)

print(f"After extinction: {len(remaining_compounds)} compounds remain")
```

### Contraction + Re-expansion (Fast Path)

For large-scale robustness analyses, use the fused method to avoid ID conversion overhead.
A single `extinctReactionSets` argument is used for both phases — no scope-intersection
needed. Contraction only operates on reactions in `reactionScope`, so extra reactions
outside scope have no effect. For re-expansion the full set is the correct mask.

```python
# 1000 random extinction sets (full, un-intersected — that's fine)
import random
all_rxns = list(kegg.rid_to_idx.keys())
extinct_sets = [random.sample(all_rxns, 500) for _ in range(1000)]

# One call: contract + re-expand, returns raw numpy arrays
result = kegg.contract_and_reexpand(seeds, reaction_scope, compound_scope, extinct_sets)

# Extract counts directly — no Python ID lookups
cpd_counts_contracted = result['x_contracted'].sum(axis=1)   # (1000,)
cpd_counts_reexpanded = result['x_reexpanded'].sum(axis=1)   # (1000,)

# Recover ID lists for specific rows when needed
o2_cpds = kegg._x_to_compounds(result['x_reexpanded'][0])
o2_rxns = kegg._y_to_reactions(result['y_reexpanded'][0])

# Null presence matrix (how many replicates include each compound)
null_cpd_presence = result['x_reexpanded'][1:].sum(axis=0)   # (n_compounds,)
```

### Trace Algorithm

Track compound and reaction addition at each iteration:

```python
compound_dict, reaction_dict = kegg.expand(seeds, algorithm='trace')

# compound_dict and reaction_dict map IDs to iteration numbers
print(f"C00001 appeared at iteration {compound_dict['C00001']}")
print(f"R00001 appeared at iteration {reaction_dict['R00001']}")
```

## API Overview

### Core Methods

#### `expand(seedSet, algorithm='naive')`
Run network expansion from seed compounds.
- **Args**:
  - `seedSet`: List of compound IDs to start from
  - `algorithm`: Either `'naive'` (fast) or `'trace'` (tracks iterations)
- **Returns**: `(compound_scope, reaction_scope)` or `(compound_dict, reaction_dict)` for trace

#### `contract(seedSet, reactionScope, compoundScope, extinctReactions)`
Contract network after removing reactions.
- **Args**:
  - `seedSet`: Seed compounds to preserve during contraction. These are compounds that entered the network as initial conditions and aren't produced by any reaction — they remain available regardless of which reactions go extinct.
  - `reactionScope`: Full reaction scope before contraction
  - `compoundScope`: Full compound scope before contraction
  - `extinctReactions`: List of reactions to remove
- **Returns**: `(remaining_compounds, remaining_reactions)`

#### `run_expansions_batch(seedSets, maskedReactionSets=None, algorithm='naive')`
Most general batch expansion — supports varying seeds and/or varying masks, parallelized via Rust/Rayon when available. Inputs are broadcast automatically: a single seed set is broadcast across all masks, a single mask across all seeds, or paired N-to-N.
- **Args**:
  - `seedSets`: A single seed set or list of seed sets
  - `maskedReactionSets`: Optional single reaction set or list of reaction sets to remove
  - `algorithm`: `'naive'` or `'trace'`
- **Returns**: `(compound_scopes, reaction_scopes)` - lists of results
- **Note**: `run_expansions()` and `run_expansions_reactionMasks()` both delegate to this method.

#### `run_expansions_reactionMasks_parallel(seedSet, maskedReactionSets)`
Run multiple masked expansions in parallel.
- **Args**:
  - `seedSet`: Initial seed compounds
  - `maskedReactionSets`: List of reaction sets to remove (one per mask)
- **Returns**: `(compound_scopes, reaction_scopes)` - lists of results

#### `run_contractions(seedSet, reactionScope, compoundScope, extinctReactionSets)`
Run multiple contractions in parallel.
- **Args**:
  - `seedSet`: Seed compounds to preserve (same across all batch rows)
  - `reactionScope`: Full reaction scope before contraction
  - `compoundScope`: Full compound scope before contraction
  - `extinctReactionSets`: List of extinction sets (each is a list of reaction IDs to remove)
- **Returns**: `(compound_scopes, reaction_scopes)` - lists of results

#### `contract_and_reexpand(seedSet, reactionScope, compoundScope, extinctReactionSets)`
Fused contraction + re-expansion in a single call, without intermediate ID conversion. Contracts the network under each extinction set, then re-expands from the surviving compounds using the same masks. All marshalling stays in numpy — ideal for large batches where you only need counts.
- **Args**:
  - `seedSet`: Seed compounds to preserve during contraction
  - `reactionScope`: Full reaction scope before contraction
  - `compoundScope`: Full compound scope before contraction
  - `extinctReactionSets`: List of extinction sets (used for both contraction and re-expansion masks). No scope-intersection needed — contraction only operates on reactions in `reactionScope`, so extra reactions outside scope have no effect. For re-expansion the full set is wanted as the mask.
- **Returns**: dict with `'x_contracted'`, `'y_contracted'`, `'x_reexpanded'`, `'y_reexpanded'` — all `(N × n_compounds/n_reactions)` uint8 numpy arrays
- **Recovering IDs**: Call `_x_to_compounds(result['x_reexpanded'][i])` or `_y_to_reactions(result['y_reexpanded'][i])` on any row you need as a list of IDs. For counts or presence matrices, operate directly on the arrays (e.g., `.sum(axis=1)`).
- **Note**: Requires the Rust backend (`netexprs`).

#### `initialize_reaction_matrix(reaction_id_sets)`
Batch version of `initialize_reaction_vector()`. Converts a list of reaction ID sets to a 2D binary matrix in one call.
- **Args**: `reaction_id_sets` — list of lists/sets of reaction IDs
- **Returns**: `(N × n_reactions)` uint8 numpy array

#### `initialize_metabolite_matrix(compound_id_sets)`
Batch version of `initialize_metabolite_vector()`. Converts a list of compound ID sets to a 2D binary matrix in one call.
- **Args**: `compound_id_sets` — list of lists/sets of compound IDs
- **Returns**: `(N × n_compounds)` uint8 numpy array

## Backend Selection

The package automatically selects the best available backend:

| Condition | Backend Used |
|-----------|--------------|
| `netexprs` installed + naive/trace algorithm | **Rust** (fast) |
| `netexprs` not installed | Python |
| `cr` or `step` algorithm requested | Python |

When Rust is active, you'll see a one-time `🦀 netexprs` message in your output on the first Rust-accelerated call.

All batch operations (multi-seed, masked, trace, contraction) use Rust with Rayon parallelism when available.

## Data Files

The package includes curated datasets in `networkExpansionPy/assets/`:

### Metabolic Networks
- `metabolism.v8.01May2023.pkl`: KEGG-based global metabolic network
- Pre-processed networks with thermodynamic constraints

### Compound Sets
- `seeds.Goldford2022.csv`: Geochemically plausible seed compounds
- `seeds.*.csv`: Various seed compound collections

## Testing

### Unit Tests
```bash
python -m pytest tests/
```

## Benchmarking

Compare Python vs Rust performance:

```bash
python benchmark/benchmark_rust_vs_python.py
```

Results are saved/appended to `benchmark_results.json`.

## Advanced Usage

### Custom Network Loading

```python
import pandas as pd

# Load custom network from DataFrame
# Required columns: 'cid' (compound), 'rn' (reaction), 's' (stoichiometry)
network_df = pd.read_csv('my_network.csv')

kegg = ne.GlobalMetabolicNetwork()
kegg.network = network_df
kegg.pruneInconsistentReactions()
```

### Working with Matrices Directly

```python
# Get stoichiometric matrices
R, P = kegg.create_RP_from_irreversible_network()
# R: reactant matrix (compounds × reactions)
# P: product matrix (compounds × reactions)

# Run low-level expansion
from scipy.sparse import csr_matrix
x0 = kegg.initialize_metabolite_vector(seeds)
b = np.sum(R, axis=0)

x_final, y_final = ne.netExp(csr_matrix(R), csr_matrix(P), 
                               csr_matrix(x0).T, csr_matrix(b).T)
```