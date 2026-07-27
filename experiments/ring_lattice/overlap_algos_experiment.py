"""
overlap_algos_experiment — compare hypercommon vs 6 baseline CD algos across overlap%.

Fixed config: n=2000, z=16, ring_size=100  (20 disjoint rings of 100).

For each overlap in {0, 2, 4, 6, 8, 10}%:
  For run_i in 1..10  (each run is a FULL fresh build):
    1. Generate G = ring_lattice + apply_overlap with fresh RNG draws (each run
       gets a different random merge pattern, not just a different rewiring trajectory).
    2. Build ground truth from the run-specific merge.
    3. Walk p = 0..1 in steps of 0.01, rewiring k_step edges at each step.
       At every p, evaluate ALL 7 algos on the same graph snapshot, in parallel.
    4. For HYPERCOMMON specifically: at every p, sweep a t-grid and record the
       BEST omega and its argmax t. The grid is full-coarse at p=0 (24 values
       step 0.04), then adaptive ±0.04 step 0.01 around the previous step's
       t_argmax for p>0 (9 values).
       This means hypercommon is evaluated at its envelope-best threshold per
       graph snapshot, not at one fixed t.

Outputs (under results/overlap_algos/run_<ts>/):
  <algo>/progress.csv                streaming append, crash-safe
  <algo>/full_sweep.xlsx             sheets o0..o10, cols p, run_1..run_10, omega_avg, omega_std
  <algo>/avg.xlsx                    one sheet 'avg', cols overlap, p, omega_avg, omega_std, n_runs
  hypercommon/t_argmax.xlsx          extra: per-(overlap, p) the t that won, per run

Notes:
  - No per-algo timeout. Hypercommon may run for a long time at high p; baselines
    are expected to finish quickly. If a baseline ever truly hangs, the experiment
    hangs — we've already dropped LFM and replaced infomap with walktrap to avoid
    the known hang/crash cases.
  - Excels are written after each overlap completes, so partial results are usable.

Run from venv:
  ./.venv/Scripts/python.exe -m experiments.ring_lattice.overlap_algos_experiment
"""

from __future__ import annotations

import os
import csv
import math
import random
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from fractions import Fraction
from statistics import mean, pstdev

import pandas as pd

# Pre-import cdlib in the main thread so its global state is initialised
# before any worker processes spawn.
warnings.filterwarnings("ignore")
from cdlib import algorithms as _cdlib_preload  # noqa: F401

from generators.ring_lattice import ring_lattice
from utils.rewiring import rewire_step
from metrics.omega import omega_index, build_pair_counts
from predicates.jaccard import closed_neighborhood_jaccard_predicate
from hypercommon.algorithm import get_communities


# =====================================================================
# Config
# =====================================================================

N         = 2000
Z         = 16
RING_SIZE = 100

OVERLAPS  = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]

P_STEP    = 0.01
AVG_FINAL = 10

# Hypercommon per-step t-grid policy.
T_INIT_STEP      = 0.04   # coarse grid step at p=0  (=> grid 0.04, 0.08, ..., 0.96)
T_ADAPTIVE_HALF  = 0.04   # +- window around prev t_argmax for p>0
T_ADAPTIVE_STEP  = 0.01   # granularity within the adaptive window

ALGOS = [
    "hypercommon",
    "leiden",
    "label_propagation",
    "walktrap",
    "slpa",
    "demon",
    "angel",
]

N_WORKERS = len(ALGOS)


# =====================================================================
# Graph build helpers (overlap merging + ground truth)
# =====================================================================

def apply_overlap(G, n: int, ring_size: int, overlap: float, rng: random.Random) -> dict[int, int]:
    """Merge floor(overlap*n) inter-ring node pairs by absorbing v into u.
    Returns {absorbed_node: surviving_node}.
    """
    k = math.floor(overlap * n)
    if k == 0:
        return {}

    node_ring = [node // ring_size for node in range(n)]
    merged: dict[int, int] = {}
    used: set[int] = set()

    attempts = 0
    max_attempts = k * 20
    while len(merged) < k and attempts < max_attempts:
        attempts += 1
        u = rng.randrange(n)
        v = rng.randrange(n)
        if u == v or node_ring[u] == node_ring[v]:
            continue
        if u in used or v in used:
            continue
        if not G.has_node(u) or not G.has_node(v):
            continue

        for neighbor in list(G.neighbors(v)):
            if neighbor != u:
                G.add_edge(u, neighbor)
        G.remove_node(v)
        used.add(v); used.add(u)
        merged[v] = u

    return merged


def ring_ground_truth_with_overlap(n: int, ring_size: int, merged: dict[int, int]) -> list[set[int]]:
    rings = n // ring_size
    communities = [set(range(r * ring_size, (r + 1) * ring_size)) for r in range(rings)]
    for v, u in merged.items():
        v_ring = v // ring_size
        communities[v_ring].discard(v)
        communities[v_ring].add(u)
    return communities


def ring_lattice_edge_count(n: int, z: int, ring_size: int) -> int:
    rings = n // ring_size
    return rings * (ring_size * z // 2)


def validate_rewiring_plan(M0: int, p_step: float):
    p = Fraction(str(p_step)).limit_denominator()
    if p.numerator != 1:
        raise ValueError(f"p_step must be 1/steps. Got {p_step}")
    steps = p.denominator
    if M0 % steps != 0:
        raise ValueError(f"M0={M0} not divisible by steps={steps}")
    return steps


def build_run_graph(overlap: float, rng: random.Random):
    """Fresh ring lattice + fresh random overlap merge. Returns (G, truth, gt_pc, n_actual, total_pairs, M_actual)."""
    rings = N // RING_SIZE
    G = ring_lattice(n=N, z=Z, rings=rings)
    merged = apply_overlap(G, n=N, ring_size=RING_SIZE, overlap=overlap, rng=rng)
    truth = ring_ground_truth_with_overlap(N, RING_SIZE, merged)
    gt_pc = build_pair_counts(truth)
    n_actual = G.number_of_nodes()
    total_pairs = n_actual * (n_actual - 1) // 2
    M_actual = G.number_of_edges()
    return G, truth, gt_pc, n_actual, total_pairs, M_actual


# =====================================================================
# Hypercommon t-grid policy
# =====================================================================

def make_t_grid(s: int, prev_t_best):
    """Coarse grid at p=0 (s=0); adaptive ±delta around prev_t_best for p>0."""
    if s == 0 or prev_t_best is None:
        # 0.04, 0.08, ..., 0.96
        n_pts = int(round((1.0 - T_INIT_STEP) / T_INIT_STEP)) + 1
        return [round((i + 1) * T_INIT_STEP, 4) for i in range(n_pts) if (i + 1) * T_INIT_STEP < 1.0]
    # local search
    half_n = int(round(T_ADAPTIVE_HALF / T_ADAPTIVE_STEP))
    grid = []
    for d in range(-half_n, half_n + 1):
        t = round(prev_t_best + d * T_ADAPTIVE_STEP, 4)
        if 0.01 <= t <= 0.99:
            grid.append(t)
    return grid


# =====================================================================
# Worker
# =====================================================================

def _run_algo_worker(algo_name, params, edges, nodes, gt_pair_counts, total_pairs):
    """Run one algo on one graph snapshot, return (algo_name, omega, t_argmax_or_None).

    For hypercommon: params = {'t_grid': [t1, t2, ...]}. Returns max omega and its t.
    For others: params = {}.
    """
    warnings.filterwarnings("ignore")

    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if algo_name == "hypercommon":
        from predicates.jaccard import closed_neighborhood_jaccard_predicate as _pred
        from hypercommon.algorithm import get_communities as _gc
        from metrics.omega import omega_index as _om, build_pair_counts as _bpc
        t_grid = params["t_grid"]
        best_omega = float("-inf")
        best_t = None
        for t in t_grid:
            try:
                pred = _gc(G, _pred(t))
                om = _om(gt_pair_counts, _bpc(pred), total_pairs)
            except Exception:
                continue
            if om > best_omega:
                best_omega = om
                best_t = t
        if best_t is None:
            return algo_name, float("nan"), None
        return algo_name, float(best_omega), float(best_t)

    try:
        if algo_name == "leiden":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.leiden(G).communities]
        elif algo_name == "label_propagation":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.label_propagation(G).communities]
        elif algo_name == "walktrap":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.walktrap(G).communities]
        elif algo_name == "slpa":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.slpa(G).communities]
        elif algo_name == "demon":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.demon(G, epsilon=0.25).communities]
        elif algo_name == "angel":
            from cdlib import algorithms as A
            pred = [set(c) for c in A.angel(G, threshold=0.25).communities]
        else:
            raise ValueError(f"unknown algo: {algo_name}")

        from metrics.omega import omega_index as _om, build_pair_counts as _bpc
        score = _om(gt_pair_counts, _bpc(pred), total_pairs)
        return algo_name, float(score), None
    except Exception:
        return algo_name, float("nan"), None


# =====================================================================
# Sweep one overlap — produces all (run, p, algo) -> (omega, t_argmax_or_None)
# =====================================================================

def sweep_one_overlap(
    overlap: float,
    rng: random.Random,
    pool: ProcessPoolExecutor,
    progress_writers: dict,
    progress_fps: dict,
    steps: int,
):
    """Returns results: dict[(run_i, step_idx, algo)] -> (omega, t_argmax_or_None)."""
    results: dict[tuple[int, int, str], tuple[float, float | None]] = {}

    for run_i in range(1, AVG_FINAL + 1):
        t_run0 = time.perf_counter()

        # ---- FULL RESET: regenerate graph + truth + edge_stack for this run ----
        G, truth, gt_pc, n_actual, total_pairs, M_actual = build_run_graph(overlap, rng)
        k_step = M_actual // steps
        edge_stack = list(G.edges())
        rng.shuffle(edge_stack)

        print(f"      [run {run_i}/{AVG_FINAL}] n_actual={n_actual} edges={M_actual} k_step={k_step}")

        prev_t_best = None  # reset per run

        for s in range(steps + 1):
            edges = list(G.edges())
            nodes = list(G.nodes())

            # Build t-grid for hypercommon (full at p=0, adaptive at p>0)
            t_grid = make_t_grid(s, prev_t_best)

            # Submit all algos in parallel
            futures = {}
            for algo in ALGOS:
                if algo == "hypercommon":
                    params = {"t_grid": t_grid}
                else:
                    params = {}
                futures[algo] = pool.submit(
                    _run_algo_worker,
                    algo, params, edges, nodes, gt_pc, total_pairs,
                )

            for algo, fut in futures.items():
                algo_name, omega, t_arg = fut.result()
                results[(run_i, s, algo_name)] = (omega, t_arg)

                row = {
                    "overlap": overlap, "run": run_i, "step": s,
                    "p": round(s / steps, 6),
                    "omega": omega,
                }
                if algo == "hypercommon":
                    row["t_argmax"] = t_arg
                    if t_arg is not None and not (isinstance(t_arg, float) and math.isnan(t_arg)):
                        prev_t_best = t_arg
                progress_writers[algo].writerow(row)

            for fp in progress_fps.values():
                fp.flush()

            if s < steps:
                rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

        dt_run = time.perf_counter() - t_run0
        print(f"      [run {run_i}/{AVG_FINAL}] done dt={dt_run:.1f}s")

    return results


# =====================================================================
# Excel writers
# =====================================================================

def overlap_label(overlap: float) -> str:
    pct = round(overlap * 100)
    return f"o{pct}"


def _omega(value):
    """Extract scalar omega from results value (which may be a (omega, t) tuple)."""
    if isinstance(value, tuple):
        return value[0]
    return value


def _t_argmax(value):
    if isinstance(value, tuple):
        return value[1]
    return None


def write_full_sweep_excels(
    root_dir: str,
    all_results: dict,
    overlaps: list[float],
    steps: int,
    avg_final: int,
):
    """Per-algo full_sweep.xlsx — one sheet per overlap, cols p, run_1..run_N, omega_avg, omega_std."""
    p_values = [round(s / steps, 6) for s in range(steps + 1)]

    for algo in ALGOS:
        algo_dir = os.path.join(root_dir, algo)
        os.makedirs(algo_dir, exist_ok=True)
        path = os.path.join(algo_dir, "full_sweep.xlsx")

        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for ov in overlaps:
                rows = []
                for s in range(steps + 1):
                    row: dict = {"p": p_values[s]}
                    run_vals = []
                    for run_i in range(1, avg_final + 1):
                        v = all_results.get((ov, run_i, s, algo))
                        omega = _omega(v) if v is not None else float("nan")
                        row[f"run_{run_i}"] = omega
                        if not (isinstance(omega, float) and math.isnan(omega)):
                            run_vals.append(omega)
                    if run_vals:
                        row["omega_avg"] = mean(run_vals)
                        row["omega_std"] = pstdev(run_vals) if len(run_vals) > 1 else 0.0
                    else:
                        row["omega_avg"] = float("nan")
                        row["omega_std"] = float("nan")
                    rows.append(row)
                df = pd.DataFrame(rows)
                df.to_excel(writer, sheet_name=overlap_label(ov), index=False)
        print(f"  [excel] {path}")


def write_avg_excels(
    root_dir: str,
    all_results: dict,
    overlaps: list[float],
    steps: int,
    avg_final: int,
):
    """Per-algo avg.xlsx — single sheet, cols overlap, p, omega_avg, omega_std, n_runs."""
    p_values = [round(s / steps, 6) for s in range(steps + 1)]

    for algo in ALGOS:
        algo_dir = os.path.join(root_dir, algo)
        os.makedirs(algo_dir, exist_ok=True)
        path = os.path.join(algo_dir, "avg.xlsx")

        rows = []
        for ov in overlaps:
            for s in range(steps + 1):
                run_vals = []
                for run_i in range(1, avg_final + 1):
                    v = all_results.get((ov, run_i, s, algo))
                    omega = _omega(v) if v is not None else float("nan")
                    if not (isinstance(omega, float) and math.isnan(omega)):
                        run_vals.append(omega)
                row = {
                    "overlap": ov,
                    "p": p_values[s],
                    "omega_avg": mean(run_vals) if run_vals else float("nan"),
                    "omega_std": (pstdev(run_vals) if len(run_vals) > 1 else 0.0) if run_vals else float("nan"),
                    "n_runs": len(run_vals),
                }
                rows.append(row)

        df = pd.DataFrame(rows)
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="avg", index=False)
        print(f"  [excel] {path}")


def write_hypercommon_t_argmax_excel(
    root_dir: str,
    all_results: dict,
    overlaps: list[float],
    steps: int,
    avg_final: int,
):
    """hypercommon/t_argmax.xlsx — sheets per overlap, cols p, t_run_1..t_run_N, t_avg."""
    p_values = [round(s / steps, 6) for s in range(steps + 1)]
    algo_dir = os.path.join(root_dir, "hypercommon")
    os.makedirs(algo_dir, exist_ok=True)
    path = os.path.join(algo_dir, "t_argmax.xlsx")

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for ov in overlaps:
            rows = []
            for s in range(steps + 1):
                row: dict = {"p": p_values[s]}
                vals = []
                for run_i in range(1, avg_final + 1):
                    v = all_results.get((ov, run_i, s, "hypercommon"))
                    t = _t_argmax(v)
                    row[f"t_run_{run_i}"] = t if t is not None else float("nan")
                    if t is not None and not (isinstance(t, float) and math.isnan(t)):
                        vals.append(t)
                row["t_avg"] = mean(vals) if vals else float("nan")
                rows.append(row)
            df = pd.DataFrame(rows)
            df.to_excel(writer, sheet_name=overlap_label(ov), index=False)
    print(f"  [excel] {path}")


# =====================================================================
# Main
# =====================================================================

def run_experiment(
    out_dir: str = "results/overlap_algos",
    seed: int = 42,
    overlaps: list[float] | None = None,
):
    if overlaps is None:
        overlaps = OVERLAPS

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    root_dir = os.path.join(out_dir, f"run_{ts}")
    os.makedirs(root_dir, exist_ok=True)

    M0 = ring_lattice_edge_count(N, Z, RING_SIZE)
    steps = validate_rewiring_plan(M0, P_STEP)

    # ---- Per-algo streaming progress.csv ----
    progress_fps: dict[str, object] = {}
    progress_writers: dict[str, csv.DictWriter] = {}
    for algo in ALGOS:
        algo_dir = os.path.join(root_dir, algo)
        os.makedirs(algo_dir, exist_ok=True)
        fp = open(os.path.join(algo_dir, "progress.csv"), "w", newline="")
        fields = ["overlap", "run", "step", "p", "omega"]
        if algo == "hypercommon":
            fields.append("t_argmax")
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        progress_fps[algo] = fp
        progress_writers[algo] = writer

    rng = random.Random(seed)
    all_results: dict = {}

    print(f"Output: {root_dir}")
    print(f"Algos: {ALGOS}")
    print(f"Overlaps: {overlaps}")
    print(f"Config: n={N}, z={Z}, ring_size={RING_SIZE}, p_step={P_STEP}, avg={AVG_FINAL}")
    print(f"Hypercommon t-grid: full coarse step {T_INIT_STEP} at p=0; adaptive +-{T_ADAPTIVE_HALF} step {T_ADAPTIVE_STEP} for p>0")
    print()

    t_global0 = time.perf_counter()

    try:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            for ov_idx, overlap in enumerate(overlaps, start=1):
                print(f"\n[OVERLAP {ov_idx}/{len(overlaps)}] overlap={overlap:.2f}")
                t_ov0 = time.perf_counter()

                results = sweep_one_overlap(
                    overlap=overlap,
                    rng=rng,
                    pool=pool,
                    progress_writers=progress_writers,
                    progress_fps=progress_fps,
                    steps=steps,
                )

                for (run_i, s, algo), v in results.items():
                    all_results[(overlap, run_i, s, algo)] = v

                # Intermediate Excel write so partial results are usable on crash
                print(f"    [excel] writing intermediate excels...")
                write_full_sweep_excels(root_dir, all_results, overlaps[:ov_idx], steps, AVG_FINAL)
                write_avg_excels(root_dir, all_results, overlaps[:ov_idx], steps, AVG_FINAL)
                write_hypercommon_t_argmax_excel(root_dir, all_results, overlaps[:ov_idx], steps, AVG_FINAL)

                dt_ov = time.perf_counter() - t_ov0
                elapsed = time.perf_counter() - t_global0
                remaining = (elapsed / ov_idx) * (len(overlaps) - ov_idx)
                print(f"[OVERLAP {ov_idx}/{len(overlaps)}] done dt={dt_ov:.1f}s  elapsed={elapsed:.1f}s  est_remaining={remaining:.1f}s")
    finally:
        for fp in progress_fps.values():
            fp.close()

    print(f"\n[ALL DONE] total dt={time.perf_counter()-t_global0:.1f}s")
    print(f"Results in: {root_dir}")


if __name__ == "__main__":
    run_experiment()
