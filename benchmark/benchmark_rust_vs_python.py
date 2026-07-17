import warnings
from scipy.sparse import SparseEfficiencyWarning
warnings.filterwarnings('ignore', category=SparseEfficiencyWarning)

import numpy as np
from scipy.sparse import csr_matrix
import networkExpansionPy.lib as ne
import netexprs
import pandas as pd
import time
import json
import os
from datetime import datetime
import multiprocessing
from random import sample, seed as random_seed

# Get the directory containing the current file
current_dir = os.path.dirname(os.path.abspath(__file__))

RESULTS_FILE = os.path.join(current_dir, "benchmark_results.json")

def load_results():
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE, "r") as f:
            return json.load(f)
    return []

def save_results(results):
    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=2)

def benchmark(fn, n_runs=10, warmup=2):
    for _ in range(warmup):
        fn()
    
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    
    return {
        "mean": np.mean(times),
        "std": np.std(times),
        "min": np.min(times),
        "max": np.max(times),
        "n_runs": n_runs,
    }

def run_benchmarks():
    random_seed(42)

    # Load data
    print("Loading data...")
    kegg = ne.load_metabolism("metabolism.v8.01May2023.pkl")
    seedSet = list(set(ne.load_compounds('seeds.Goldford2022.csv')["ID"]))
    
    kegg._ensure_dicts()
    kegg._ensure_rust_ready()
    ra = kegg._rust_arrays

    R, P = kegg.create_RP_from_irreversible_network()
    b_dense = sum(R)
    R = csr_matrix(R)
    P = csr_matrix(P)
    b_sparse = csr_matrix(b_dense).transpose()

    x0 = kegg.initialize_metabolite_vector(seedSet)
    x0_sparse = csr_matrix(x0).transpose()

    # Arrays for Rust calls
    n_reactions = ra["n_reactions"]
    n_compounds = ra["n_compounds"]
    x0_arr = x0.astype(np.uint8)
    x0_2d = x0_arr.reshape(1, -1)

    all_rxns = list(kegg.rid_to_idx.keys())
    network_rxns = list(set(r[0] for r in all_rxns))

    results = {
        "timestamp": datetime.now().isoformat(),
        "n_compounds": n_compounds,
        "n_reactions": n_reactions,
        "n_cores": multiprocessing.cpu_count(),
        "benchmarks": {}
    }
    
    # ── 1. Single Expansion (no mask) ────────────────────────────────────
    print("\nBenchmarking single expansion...")

    results["benchmarks"]["expand_python"] = benchmark(
        lambda: ne.netExp(R, P, x0_sparse, b_sparse),
        n_runs=10,
    )
    results["benchmarks"]["expand_rust"] = benchmark(
        lambda: netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x0_2d.copy(), ra["b"], None,
        ),
        n_runs=10,
    )

    # ── 2. Single Masked Expansion ───────────────────────────────────────
    print("Benchmarking single masked expansion...")

    np.random.seed(42)
    mask_1d = (np.random.random(n_reactions) > 0.1).astype(np.uint8)
    mask_2d = mask_1d.reshape(1, -1)

    results["benchmarks"]["expand_masked_python"] = benchmark(
        lambda: ne.netExp(
            R * csr_matrix(np.diag(mask_1d)),
            P * csr_matrix(np.diag(mask_1d)),
            x0_sparse, b_sparse,
        ),
        n_runs=10,
    )
    results["benchmarks"]["expand_masked_rust"] = benchmark(
        lambda: netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x0_2d.copy(), ra["b"], mask_2d,
        ),
        n_runs=10,
    )

    # ── 3. Batch Masked Expansion (n=10) ─────────────────────────────────
    n_masks = 10
    print(f"Benchmarking batch masked expansion (n={n_masks})...")

    np.random.seed(123)
    masks = (np.random.random((n_masks, n_reactions)) > 0.1).astype(np.uint8)
    x_tiled = np.tile(x0_arr, (n_masks, 1))

    # Python: convert masks to reaction ID lists for the high-level wrapper
    masked_reaction_sets = []
    for i in range(n_masks):
        removed_rxns = [kegg.idx_to_rid[j] for j in range(n_reactions) if masks[i, j] == 0]
        masked_reaction_sets.append(removed_rxns)

    results["benchmarks"]["expand_masked_batch_10_python"] = benchmark(
        lambda: kegg._run_expansions_reactionMasks_python(seedSet, masked_reaction_sets),
        n_runs=3,
    )
    results["benchmarks"]["expand_masked_batch_10_rust"] = benchmark(
        lambda: netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x_tiled.copy(), ra["b"], masks,
        ),
        n_runs=3,
    )

    # ── 4. Batch Multi-Seed Expansion (n=10) ─────────────────────────────
    n_seeds = 10
    print(f"Benchmarking batch multi-seed expansion (n={n_seeds})...")

    random_seed(99)
    seed_sets = []
    x_seed_batch = np.zeros((n_seeds, n_compounds), dtype=np.uint8)
    for i in range(n_seeds):
        s = sample(seedSet, max(1, len(seedSet) - i * 5))
        seed_sets.append(s)
        x_seed_batch[i] = kegg.initialize_metabolite_vector(s).astype(np.uint8)

    results["benchmarks"]["expand_multiseed_10_python"] = benchmark(
        lambda: kegg._run_expansions_python(seed_sets, "naive"),
        n_runs=3,
    )
    results["benchmarks"]["expand_multiseed_10_rust"] = benchmark(
        lambda: netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x_seed_batch.copy(), ra["b"], None,
        ),
        n_runs=3,
    )

    # ── 5. Single Trace Expansion ────────────────────────────────────────
    print("Benchmarking single trace expansion...")

    results["benchmarks"]["trace_python"] = benchmark(
        lambda: ne.netExp_trace(R, P, x0_sparse, b_sparse),
        n_runs=10,
    )
    results["benchmarks"]["trace_rust"] = benchmark(
        lambda: netexprs.expand_trace_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x0_2d.copy(), ra["b"], None,
        ),
        n_runs=10,
    )

    # ── 6. Batch Trace Expansion (n=10) ──────────────────────────────────
    print(f"Benchmarking batch trace expansion (n={n_seeds})...")

    results["benchmarks"]["trace_batch_10_python"] = benchmark(
        lambda: kegg._run_expansions_python(seed_sets, "trace"),
        n_runs=3,
    )
    results["benchmarks"]["trace_batch_10_rust"] = benchmark(
        lambda: netexprs.expand_trace_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x_seed_batch.copy(), ra["b"], None,
        ),
        n_runs=3,
    )

    # ── Setup for Contraction ────────────────────────────────────────────
    print("\nRunning initial expansion for contraction benchmarks...")
    x_expanded_sparse, y_expanded_sparse = ne.netExp(R, P, x0_sparse, b_sparse)

    x_expanded_arr = np.asarray(x_expanded_sparse.toarray()).ravel().astype(np.uint8)
    y_expanded_arr = np.asarray(y_expanded_sparse.toarray()).ravel().astype(np.uint8)

    cidx = np.nonzero(x_expanded_sparse.toarray().T[0])[0]
    compoundScope = [kegg.idx_to_cid[i] for i in cidx]
    ridx = np.nonzero(y_expanded_sparse.toarray().T[0])[0]
    reactionScope = [kegg.idx_to_rid[i] for i in ridx]

    # ── 7. Single Contraction ────────────────────────────────────────────
    print("Benchmarking single contraction...")

    np.random.seed(55)
    y_extinct_arr = (np.random.random(n_reactions) < 0.05).astype(np.uint8)
    y_extinct_sparse = csr_matrix(y_extinct_arr).T

    results["benchmarks"]["contract_python"] = benchmark(
        lambda: ne.netContract(R, P, b_sparse, x_expanded_sparse, y_expanded_sparse, y_extinct_sparse),
        n_runs=10,
    )
    results["benchmarks"]["contract_rust"] = benchmark(
        lambda: netexprs.contract_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            x_expanded_arr.reshape(1, -1),
            y_expanded_arr.reshape(1, -1),
            y_extinct_arr.reshape(1, -1),
        ),
        n_runs=10,
    )

    # ── 8. Batch Contraction (n=10) ──────────────────────────────────────
    n_batches = 10
    print(f"Benchmarking batch contraction (n={n_batches})...")

    np.random.seed(999)
    y_extinct_batch = (np.random.random((n_batches, n_reactions)) < 0.05).astype(np.uint8)

    extinct_reaction_sets = []
    for i in range(n_batches):
        extinct_indices = np.where(y_extinct_batch[i] == 1)[0]
        extinct_ids = [kegg.idx_to_rid[idx] for idx in extinct_indices]
        extinct_reaction_sets.append(extinct_ids)

    results["benchmarks"]["contract_batch_10_python"] = benchmark(
        lambda: kegg._run_contractions_python(reactionScope, compoundScope, extinct_reaction_sets),
        n_runs=3,
    )
    results["benchmarks"]["contract_batch_10_rust"] = benchmark(
        lambda: netexprs.contract_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], n_reactions,
            ra["p_data"], ra["p_indices"], ra["p_indptr"], n_compounds,
            np.tile(x_expanded_arr, (n_batches, 1)),
            np.tile(y_expanded_arr, (n_batches, 1)),
            y_extinct_batch,
        ),
        n_runs=3,
    )

    # ── Print Summary ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Cores: {results['n_cores']}")
    print(f"Compounds: {results['n_compounds']}, Reactions: {results['n_reactions']}")
    print()

    for name, stats in results["benchmarks"].items():
        print(f"{name:<40}: {stats['mean']*1000:>9.2f} ms (±{stats['std']*1000:.2f} ms)")

    print("\n" + "-" * 60)
    print("SPEEDUPS (Python / Rust)")
    print("-" * 60)

    pairs = [
        ("expand_python", "expand_rust", "Single Expansion"),
        ("expand_masked_python", "expand_masked_rust", "Single Masked Expansion"),
        ("expand_masked_batch_10_python", "expand_masked_batch_10_rust", "Batch Masked Exp (n=10)"),
        ("expand_multiseed_10_python", "expand_multiseed_10_rust", "Batch Multi-Seed Exp (n=10)"),
        ("trace_python", "trace_rust", "Single Trace"),
        ("trace_batch_10_python", "trace_batch_10_rust", "Batch Trace (n=10)"),
        ("contract_python", "contract_rust", "Single Contraction"),
        ("contract_batch_10_python", "contract_batch_10_rust", "Batch Contraction (n=10)"),
    ]

    for py, rs, label in pairs:
        if py in results["benchmarks"] and rs in results["benchmarks"]:
            speedup = results["benchmarks"][py]["mean"] / results["benchmarks"][rs]["mean"]
            print(f"{label:<35}: {speedup:>7.1f}x")

    # Save results
    all_results = load_results()
    all_results.append(results)
    save_results(all_results)
    print(f"\nResults saved to {RESULTS_FILE}")

if __name__ == "__main__":
    run_benchmarks()