"""
pcrit_tbest_omega — full omega(p, t) surface for hypercommon at one fixed config.

Config: n=2000, z=16, ring_size=100, overlap=0  (clean ring lattice, no merging).

For each of N_RUNS independent rewiring trajectories:
  Walk p = 0..1 in steps of 0.01 (101 values). At each p:
    Sweep t in {0.00, 0.01, ..., 1.00} (101 values).
    Record omega(p, t).
  Then rewire k_step edges and continue.

Per-run outputs (under results/pcrit_tbest_omega/run_<ts>/run_<i>/):
  sweep.csv            streaming, crash-safe; cols p, t, omega
  argmax.csv           cols p, t_argmax, omega_at_argmax
  grid.xlsx            pivot: rows = p, cols = t (101 columns), cells = omega
                       written after the run completes

This script does NO averaging. We save the 3 runs separately so plotting code
can choose how to aggregate (mean / median / per-run overlay / etc).

Run from venv:
  ./.venv/Scripts/python.exe -m experiments.ring_lattice.pcrit_tbest_omega
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

import pandas as pd

warnings.filterwarnings("ignore")

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
OVERLAP   = 0.0      # this experiment is fixed at zero overlap

P_STEP    = 0.01     # p in {0.00, 0.01, ..., 1.00} -> 101 values
T_STEP    = 0.01     # t in {0.00, 0.01, ..., 1.00} -> 101 values

N_RUNS    = 3
N_WORKERS = 7        # parallel workers for the t-sweep within each p


def t_grid() -> list[float]:
    """t in {0.00, 0.01, ..., 1.00}."""
    return [round(i * T_STEP, 4) for i in range(int(round(1.0 / T_STEP)) + 1)]


# =====================================================================
# Worker
# =====================================================================

def _hypercommon_omega_at_t(t, edges, nodes, gt_pair_counts, total_pairs):
    """Run hypercommon at threshold t on the given graph snapshot. Returns (t, omega)."""
    import warnings as _w; _w.filterwarnings("ignore")
    import networkx as nx
    from predicates.jaccard import closed_neighborhood_jaccard_predicate as _pred
    from hypercommon.algorithm import get_communities as _gc
    from metrics.omega import omega_index as _om, build_pair_counts as _bpc

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    try:
        pred = _gc(G, _pred(t))
        omega = _om(gt_pair_counts, _bpc(pred), total_pairs)
        return float(t), float(omega)
    except Exception:
        return float(t), float("nan")


# =====================================================================
# Per-run sweep
# =====================================================================

def run_one(run_i: int, run_dir: str, rng: random.Random, pool: ProcessPoolExecutor, steps: int) -> None:
    os.makedirs(run_dir, exist_ok=True)

    # ---- Open streaming CSVs ----
    sweep_path  = os.path.join(run_dir, "sweep.csv")
    argmax_path = os.path.join(run_dir, "argmax.csv")

    sweep_fp = open(sweep_path, "w", newline="")
    sweep_w  = csv.DictWriter(sweep_fp, fieldnames=["p", "t", "omega"])
    sweep_w.writeheader()

    argmax_fp = open(argmax_path, "w", newline="")
    argmax_w  = csv.DictWriter(argmax_fp, fieldnames=["p", "t_argmax", "omega_at_argmax"])
    argmax_w.writeheader()

    # ---- Build graph (overlap=0, so just a ring lattice) ----
    rings = N // RING_SIZE
    G = ring_lattice(n=N, z=Z, rings=rings)
    truth = [set(range(r * RING_SIZE, (r + 1) * RING_SIZE)) for r in range(rings)]
    gt_pc = build_pair_counts(truth)
    total_pairs = N * (N - 1) // 2
    M_actual = G.number_of_edges()
    k_step = M_actual // steps

    edge_stack = list(G.edges())
    rng.shuffle(edge_stack)

    print(f"  [run {run_i}/{N_RUNS}] n={N} z={Z} edges={M_actual} k_step={k_step}")

    TGRID = t_grid()
    grid_rows: list[dict] = []

    t_run0 = time.perf_counter()

    try:
        for s in range(steps + 1):
            p = round(s / steps, 6)
            edges = list(G.edges())
            nodes = list(G.nodes())
            t_p0 = time.perf_counter()

            # Submit all t-evals; collect in t-order
            futures = {t: pool.submit(_hypercommon_omega_at_t, t, edges, nodes, gt_pc, total_pairs) for t in TGRID}

            row_omega: dict[float, float] = {}
            for t in TGRID:
                _, omega = futures[t].result()
                row_omega[t] = omega
                sweep_w.writerow({"p": p, "t": t, "omega": omega})
            sweep_fp.flush()

            # argmax
            valid = {t: o for t, o in row_omega.items() if not (isinstance(o, float) and math.isnan(o))}
            if valid:
                t_arg = max(valid, key=valid.get)
                argmax_w.writerow({"p": p, "t_argmax": t_arg, "omega_at_argmax": valid[t_arg]})
            else:
                argmax_w.writerow({"p": p, "t_argmax": "", "omega_at_argmax": ""})
            argmax_fp.flush()

            # grid row for this p (one column per t)
            row = {"p": p}
            for t in TGRID:
                row[f"t={t:.2f}"] = row_omega[t]
            grid_rows.append(row)

            dt_p = time.perf_counter() - t_p0
            elapsed = time.perf_counter() - t_run0
            est_remaining = (elapsed / (s + 1)) * (steps - s)
            if (s + 1) % 5 == 0 or s == 0:
                print(f"    [run {run_i}] p={p:.2f}  dt_p={dt_p:.1f}s  elapsed={elapsed:.1f}s  est_remaining={est_remaining:.1f}s")

            if s < steps:
                rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)
    finally:
        sweep_fp.close()
        argmax_fp.close()

    # ---- grid.xlsx (one sheet, rows=p, cols=t) ----
    grid_path = os.path.join(run_dir, "grid.xlsx")
    grid_df = pd.DataFrame(grid_rows)
    with pd.ExcelWriter(grid_path, engine="openpyxl") as writer:
        grid_df.to_excel(writer, sheet_name="omega", index=False)

    dt_run = time.perf_counter() - t_run0
    print(f"  [run {run_i}/{N_RUNS}] done dt={dt_run:.1f}s  -> {run_dir}")


# =====================================================================
# Main
# =====================================================================

def ring_lattice_edge_count(n: int, z: int, ring_size: int) -> int:
    rings = n // ring_size
    return rings * (ring_size * z // 2)


def validate_rewiring_plan(M0: int, p_step: float) -> int:
    p = Fraction(str(p_step)).limit_denominator()
    if p.numerator != 1:
        raise ValueError(f"p_step must be 1/steps. Got {p_step}")
    steps = p.denominator
    if M0 % steps != 0:
        raise ValueError(f"M0={M0} not divisible by steps={steps}")
    return steps


def run_experiment(out_dir: str = "results/pcrit_tbest_omega", seed: int = 42, n_runs: int = N_RUNS):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    root_dir = os.path.join(out_dir, f"run_{ts}")
    os.makedirs(root_dir, exist_ok=True)

    M0 = ring_lattice_edge_count(N, Z, RING_SIZE)
    steps = validate_rewiring_plan(M0, P_STEP)

    rng = random.Random(seed)

    t_global0 = time.perf_counter()
    print(f"Output: {root_dir}")
    print(f"Config: n={N}, z={Z}, ring_size={RING_SIZE}, overlap={OVERLAP}")
    print(f"p in [0,1] step {P_STEP} ({steps + 1} values)")
    print(f"t in [0,1] step {T_STEP} ({len(t_grid())} values)")
    print(f"runs={n_runs}, workers={N_WORKERS}")
    print()

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for run_i in range(1, n_runs + 1):
            run_dir = os.path.join(root_dir, f"run_{run_i}")
            run_one(run_i=run_i, run_dir=run_dir, rng=rng, pool=pool, steps=steps)

    print(f"\n[ALL DONE] total dt={time.perf_counter()-t_global0:.1f}s")
    print(f"Results in: {root_dir}")


if __name__ == "__main__":
    run_experiment()
