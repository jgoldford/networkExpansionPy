"""
MWE to profile contraction + re-expansion pipeline.

Run with: python debug_pipeline_timing.py

Identifies where time is spent:
  - Python-side vector construction (initialize_reaction_vector loops)
  - Rust contraction
  - Python-side vector construction for re-expansion
  - Rust re-expansion
"""

import time
import numpy as np
import pickle
import networkExpansionPy.lib as ne

# ── Setup ────────────────────────────────────────────────────────────────

print("Loading metabolism...")
kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
seeds = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))

kegg._ensure_dicts()
kegg._ensure_rust_ready()
ra = kegg._rust_arrays

print(f"  n_compounds={ra['n_compounds']}, n_reactions={ra['n_reactions']}")
print(f"  n_seeds={len(seeds)}")

# ── Get a realistic scope (full expansion) ───────────────────────────────

print("\nRunning initial expansion...")
cpd_scope, rxn_scope = kegg._expand_rust(seeds)
print(f"  Scope: {len(cpd_scope)} compounds, {len(rxn_scope)} reactions")

# rxn_scope contains (rn_id, direction) tuples
# Extract plain rn ID strings (like the notebook does)
scope_rn_strings = list(set(r[0] for r in rxn_scope))
scope_rn_tuples = list(rxn_scope)

print(f"  Unique rn strings: {len(scope_rn_strings)}")
print(f"  Directed tuples: {len(scope_rn_tuples)}")

# ── Build fake extinct sets (mimicking stratified null) ──────────────────

N_NULL = 1000
np.random.seed(42)
n_remove = len(scope_rn_strings) // 10  # ~10% removal rate

print(f"\nBuilding {N_NULL} extinct sets, each removing ~{n_remove} rn IDs...")

t0 = time.perf_counter()
extinct_sets_strings = []
for _ in range(N_NULL):
    idx = np.random.choice(len(scope_rn_strings), n_remove, replace=False)
    extinct_sets_strings.append([scope_rn_strings[i] for i in idx])
t1 = time.perf_counter()
print(f"  Building extinct ID lists: {t1-t0:.3f}s")

# Also build tuple versions for comparison
t0 = time.perf_counter()
extinct_sets_tuples = []
for es in extinct_sets_strings:
    tuples = []
    for rn in es:
        tuples.append((rn, "forward"))
        tuples.append((rn, "reverse"))
    extinct_sets_tuples.append(tuples)
t1 = time.perf_counter()
print(f"  Converting to tuple lists: {t1-t0:.3f}s")


# ── Profile: initialize_reaction_vector ──────────────────────────────────

print("\n" + "="*60)
print("PROFILING: initialize_reaction_vector")
print("="*60)

# Single call with string IDs
t0 = time.perf_counter()
for _ in range(10):
    kegg.initialize_reaction_vector(extinct_sets_strings[0])
t1 = time.perf_counter()
print(f"  Single call (strings, {len(extinct_sets_strings[0])} IDs): {(t1-t0)/10*1000:.2f}ms")

# Single call with tuple IDs
t0 = time.perf_counter()
for _ in range(10):
    kegg.initialize_reaction_vector(extinct_sets_tuples[0])
t1 = time.perf_counter()
print(f"  Single call (tuples, {len(extinct_sets_tuples[0])} IDs): {(t1-t0)/10*1000:.2f}ms")

# N_NULL calls with string IDs (this is what run_contractions does internally)
t0 = time.perf_counter()
for es in extinct_sets_strings:
    kegg.initialize_reaction_vector(es)
t1 = time.perf_counter()
print(f"  {N_NULL} calls (strings): {t1-t0:.3f}s")

# N_NULL calls with tuple IDs
t0 = time.perf_counter()
for es in extinct_sets_tuples:
    kegg.initialize_reaction_vector(es)
t1 = time.perf_counter()
print(f"  {N_NULL} calls (tuples): {t1-t0:.3f}s")

# Direct numpy construction (bypass initialize_reaction_vector entirely)
t0 = time.perf_counter()
y_extinct_batch = np.zeros((N_NULL, ra["n_reactions"]), dtype=np.uint8)
for i, es in enumerate(extinct_sets_strings):
    for rn in es:
        for direction in ("forward", "reverse"):
            key = (rn, direction)
            if key in kegg.rid_to_idx:
                y_extinct_batch[i, kegg.rid_to_idx[key]] = 1
t1 = time.perf_counter()
print(f"  {N_NULL} direct numpy builds (strings): {t1-t0:.3f}s")


# ── Profile: run_contractions (full pipeline) ────────────────────────────

print("\n" + "="*60)
print("PROFILING: run_contractions")
print("="*60)

# With string IDs
t0 = time.perf_counter()
cpd_contracted, rxn_contracted = kegg.run_contractions(
    None, scope_rn_strings, cpd_scope, extinct_sets_strings[:100]
)
t1 = time.perf_counter()
print(f"  run_contractions (100 string extinct sets): {t1-t0:.3f}s")

# Break it down: vector construction vs Rust call
t0 = time.perf_counter()
x_active = kegg.initialize_metabolite_vector(cpd_scope).astype(np.uint8)
y_active = kegg.initialize_reaction_vector(scope_rn_strings).astype(np.uint8)
t1 = time.perf_counter()
print(f"    Scope vector construction: {t1-t0:.3f}s")

t0 = time.perf_counter()
y_extinct_batch = np.zeros((100, ra["n_reactions"]), dtype=np.uint8)
for i, es in enumerate(extinct_sets_strings[:100]):
    y_extinct_batch[i] = kegg.initialize_reaction_vector(es).astype(np.uint8)
t1 = time.perf_counter()
print(f"    100 extinct vectors (strings): {t1-t0:.3f}s")

t0 = time.perf_counter()
x_batch, y_batch = kegg._call_contract_batch(
    np.tile(x_active, (100, 1)),
    np.tile(y_active, (100, 1)),
    y_extinct_batch,
)
t1 = time.perf_counter()
print(f"    Rust contract_batch (100): {t1-t0:.3f}s")

t0 = time.perf_counter()
for i in range(100):
    kegg._x_to_compounds(x_batch[i])
    kegg._y_to_reactions(y_batch[i])
t1 = time.perf_counter()
print(f"    Output conversion (100): {t1-t0:.3f}s")


# ── Profile: run_expansions_batch (full pipeline) ────────────────────────

print("\n" + "="*60)
print("PROFILING: run_expansions_batch")
print("="*60)

# Use contracted results as seeds, extinct sets as masks
post_contraction_cpds = cpd_contracted[:100]

t0 = time.perf_counter()
reexp_cpds, reexp_rxns = kegg.run_expansions_batch(
    post_contraction_cpds, extinct_sets_strings[:100]
)
t1 = time.perf_counter()
print(f"  run_expansions_batch (100 paired): {t1-t0:.3f}s")

# Break it down
t0 = time.perf_counter()
x_init_batch = np.zeros((100, ra["n_compounds"]), dtype=np.uint8)
for i, cpds in enumerate(post_contraction_cpds):
    x_init_batch[i] = kegg.initialize_metabolite_vector(cpds).astype(np.uint8)
t1 = time.perf_counter()
print(f"    100 seed vectors: {t1-t0:.3f}s")

t0 = time.perf_counter()
masks = np.ones((100, ra["n_reactions"]), dtype=np.uint8)
for i, rxns_removed in enumerate(extinct_sets_strings[:100]):
    exclude_vec = kegg.initialize_reaction_vector(rxns_removed).astype(np.uint8)
    masks[i] = 1 - exclude_vec
t1 = time.perf_counter()
print(f"    100 mask vectors (strings): {t1-t0:.3f}s")

t0 = time.perf_counter()
x_batch, y_batch = kegg._call_expand_batch(x_init_batch, masks)
t1 = time.perf_counter()
print(f"    Rust expand_batch (100): {t1-t0:.3f}s")

t0 = time.perf_counter()
for i in range(100):
    kegg._x_to_compounds(x_batch[i])
    kegg._y_to_reactions(y_batch[i])
t1 = time.perf_counter()
print(f"    Output conversion (100): {t1-t0:.3f}s")


# ── Profile: full pipeline at N_NULL scale ───────────────────────────────

print("\n" + "="*60)
print(f"PROFILING: full pipeline at N_NULL={N_NULL}")
print("="*60)

t_total_start = time.perf_counter()

# Contraction
t0 = time.perf_counter()
cpd_contracted_full, rxn_contracted_full = kegg.run_contractions(
    None, scope_rn_strings, cpd_scope, extinct_sets_strings
)
t1 = time.perf_counter()
print(f"  run_contractions ({N_NULL}): {t1-t0:.3f}s")

# Re-expansion
t0 = time.perf_counter()
reexp_cpds_full, reexp_rxns_full = kegg.run_expansions_batch(
    cpd_contracted_full, extinct_sets_strings
)
t1 = time.perf_counter()
print(f"  run_expansions_batch ({N_NULL}): {t1-t0:.3f}s")

t_total = time.perf_counter() - t_total_start
print(f"\n  TOTAL for one iteration: {t_total:.3f}s")
print(f"  Estimated for 70 iterations: {t_total*70:.0f}s ({t_total*70/60:.1f}min)")