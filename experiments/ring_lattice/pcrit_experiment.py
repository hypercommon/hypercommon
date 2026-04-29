"""
p_crit experiment — find the critical rewiring probability for community detection.

For each (overlap, n, z, ring_size) config:
  1. Coarse threshold search: step t by 0.04 from 0.04 upward, evaluate p_crit(t)
     for each t, stop as soon as we pass the peak (previous was better than current).
  2. Fine threshold search: step t by 0.01 in [t_peak - 0.03, t_peak + 0.03],
     evaluate p_crit(t) for each, pick t_best.
  3. Full p-sweep at t_best with avg_times=10, save everything.

p_crit(t) is estimated by walking p from 0->1 and stopping at first crossing of 0.6.
All data is saved to CSV after every evaluation so nothing is lost on crash.

Overlap: floor(overlap * n) inter-ring node pairs are merged (v absorbed into u).
The surviving node u appears in both communities in the ground truth.
"""

import os
import csv
import random
import math
import time
from datetime import datetime
from fractions import Fraction
from itertools import product

from generators.ring_lattice import ring_lattice
from utils.rewiring import rewire_step
from metrics.omega import omega_index, build_pair_counts
from predicates.jaccard import closed_neighborhood_jaccard_predicate
from hypercommon.algorithm import get_communities


# ---------------------------------------------------------------------------
# Parameter space
# ---------------------------------------------------------------------------

OVERLAPS    = [0.01]
NS          = [400, 800, 1200, 2000, 3000]
ZS          = [4, 8, 16, 32]
RING_SIZES  = [50, 100, 200]

P_STEP          = 0.01
AVG_SEARCH      = 5
AVG_FINAL       = 10
P_CRIT_THRESH   = 0.6

COARSE_T_STEP   = 0.04
FINE_T_STEP     = 0.01
FINE_T_WINDOW   = 0.03   # +- around coarse peak


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_valid_config(n, z, ring_size):
    return (
        n % ring_size == 0
        and z % 2 == 0
        and z < ring_size
    )


def all_configs():
    """Ordered: overlap -> n -> z -> ring_size."""
    configs = []
    for overlap, n, z, ring_size in product(OVERLAPS, NS, ZS, RING_SIZES):
        if is_valid_config(n, z, ring_size):
            configs.append((overlap, n, z, ring_size))
    return configs


def remaining_configs():
    """The 4 configs that did not finish in run_20260423_221250 due to crash."""
    return [
        (0.01, 3000, 16, 200),
        (0.01, 3000, 32,  50),
        (0.01, 3000, 32, 100),
        (0.01, 3000, 32, 200),
    ]


def apply_overlap(G, n: int, ring_size: int, overlap: float, rng: random.Random) -> dict[int, int]:
    """
    Merge floor(overlap * n) inter-ring node pairs by absorbing v into u.

    For each chosen pair (u, v) from different rings:
      - all edges of v are redirected to u (self-loops skipped)
      - v is removed from G

    Returns merged: {absorbed_node: surviving_node}
    """
    k = math.floor(overlap * n)
    if k == 0:
        return {}

    # ring membership for each node (only original nodes, before any removal)
    node_ring = [node // ring_size for node in range(n)]

    merged = {}      # absorbed -> surviving
    used = set()     # nodes already involved in any merge (either role)

    attempts = 0
    max_attempts = k * 20

    while len(merged) < k and attempts < max_attempts:
        attempts += 1
        u = rng.randrange(n)
        v = rng.randrange(n)

        if u == v:
            continue
        if node_ring[u] == node_ring[v]:
            continue
        if u in used or v in used:
            continue
        if not G.has_node(u) or not G.has_node(v):
            continue

        # absorb v into u
        for neighbor in list(G.neighbors(v)):
            if neighbor != u:
                G.add_edge(u, neighbor)
        G.remove_node(v)
        used.add(v)
        used.add(u)
        merged[v] = u

    return merged


def ring_ground_truth_with_overlap(n: int, ring_size: int, merged: dict[int, int]) -> list[set[int]]:
    """
    Build overlapping ground truth communities.

    Base: each ring is one community containing its original node range.
    For each (absorbed v -> surviving u):
      - remove v from its ring community (v no longer exists as a node)
      - add u to v's ring community (u now represents v, so u is in both rings)
    """
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
        raise ValueError(f"p_step must be of the form 1/steps. Got {p_step}")
    steps = p.denominator
    if M0 % steps != 0:
        raise ValueError(f"M0={M0} not divisible by steps={steps} for p_step={p_step}")
    return steps, M0 // steps


def frange(start, stop, step):
    """Inclusive float range rounded to avoid drift."""
    vals = []
    i = 0
    while True:
        v = round(start + i * step, 6)
        if v > round(stop, 6) + 1e-9:
            break
        vals.append(v)
        i += 1
    return vals


def estimate_p_crit(G_base, gt_pair_counts, total_pairs, threshold, steps, k_step, avg_times, rng):
    """
    Estimate p_crit for a given threshold by averaging avg_times independent
    forward walks. Each walk goes p=0->1 step by step and stops at first crossing
    of P_CRIT_THRESH. Returns the average p at which score dropped below threshold,
    or 1.0 if it never dropped.
    """
    pred_fn = closed_neighborhood_jaccard_predicate(threshold)
    total_p_crit = 0.0

    for _ in range(avg_times):
        G = G_base.copy()
        edge_stack = list(G.edges())
        rng.shuffle(edge_stack)

        found = False
        for s in range(steps + 1):
            pred = get_communities(G, pred_fn)
            pred_pair_count = build_pair_counts(pred)
            score = omega_index(
                ground_truth_pair_counts=gt_pair_counts,
                detected_pair_counts=pred_pair_count,
                total_pairs=total_pairs,
            )
            if score < P_CRIT_THRESH:
                total_p_crit += round(s / steps, 6)
                found = True
                break
            if s < steps:
                rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

        if not found:
            total_p_crit += 1.0

    return total_p_crit / avg_times


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def write_search_row(writer, overlap, n, z, ring_size, phase, threshold, p_crit_est):
    writer.writerow({
        "overlap": overlap, "n": n, "z": z, "ring_size": ring_size,
        "phase": phase, "threshold": threshold, "p_crit_est": p_crit_est,
    })


def write_sweep_row(writer, overlap, n, z, ring_size, t_best, p, score):
    writer.writerow({
        "overlap": overlap, "n": n, "z": z, "ring_size": ring_size,
        "t_best": t_best, "p": p, "score": score,
    })


def write_summary_row(writer, overlap, n, z, ring_size, t_best, p_crit):
    writer.writerow({
        "overlap": overlap, "n": n, "z": z, "ring_size": ring_size,
        "t_best": t_best,
        "p_crit": p_crit if p_crit is not None else "N/A",
    })


# ---------------------------------------------------------------------------
# Per-config experiment
# ---------------------------------------------------------------------------

def run_config(overlap, n, z, ring_size, rng, search_writer, sweep_writer, summary_writer):
    rings = n // ring_size

    G_base = ring_lattice(n=n, z=z, rings=rings)
    merged = apply_overlap(G_base, n=n, ring_size=ring_size, overlap=overlap, rng=rng)

    truth = ring_ground_truth_with_overlap(n=n, ring_size=ring_size, merged=merged)
    gt_pair_counts = build_pair_counts(truth)

    n_actual = G_base.number_of_nodes()
    total_pairs = n_actual * (n_actual - 1) // 2

    M0 = ring_lattice_edge_count(n, z, ring_size)
    steps, _ = validate_rewiring_plan(M0, P_STEP)
    # k_step derived from actual post-merge edge count so edge_stack never runs dry
    M_actual = G_base.number_of_edges()
    k_step = M_actual // steps

    # ---- Phase 1: coarse t search — stop after peak ----
    print(f"    [SEARCH coarse] step={COARSE_T_STEP}, stopping after peak")
    coarse_results = []
    t = COARSE_T_STEP

    while t < 1.0:
        t0 = time.perf_counter()
        p_est = estimate_p_crit(G_base, gt_pair_counts, total_pairs, t, steps, k_step, AVG_SEARCH, rng)
        coarse_results.append((t, p_est))
        write_search_row(search_writer, overlap, n, z, ring_size, "coarse", t, p_est)
        print(f"      t={t:.2f}  p_crit_est={p_est:.4f}  dt={time.perf_counter()-t0:.2f}s")

        if len(coarse_results) >= 2 and coarse_results[-1][1] < coarse_results[-2][1]:
            break

        t = round(t + COARSE_T_STEP, 6)

    best_coarse_t, best_coarse_p = max(coarse_results, key=lambda x: x[1])
    print(f"    [SEARCH coarse] peak t={best_coarse_t:.2f}  p_crit_est={best_coarse_p:.4f}")

    # ---- Phase 2: fine t search around peak ----
    fine_min = max(0.01, round(best_coarse_t - FINE_T_WINDOW, 6))
    fine_max = min(0.99, round(best_coarse_t + FINE_T_WINDOW, 6))
    print(f"    [SEARCH fine]   t in [{fine_min}, {fine_max}] step={FINE_T_STEP}")
    fine_results = []

    for t in frange(fine_min, fine_max, FINE_T_STEP):
        t0 = time.perf_counter()
        p_est = estimate_p_crit(G_base, gt_pair_counts, total_pairs, t, steps, k_step, AVG_SEARCH, rng)
        fine_results.append((t, p_est))
        write_search_row(search_writer, overlap, n, z, ring_size, "fine", t, p_est)
        print(f"      t={t:.2f}  p_crit_est={p_est:.4f}  dt={time.perf_counter()-t0:.2f}s")

    t_best, _ = max(fine_results, key=lambda x: x[1])
    print(f"    [SEARCH fine]   t_best={t_best:.2f}")

    # ---- Phase 3: full p-sweep at t_best ----
    print(f"    [SWEEP] t_best={t_best:.2f}, p 0->1 step={P_STEP}, avg={AVG_FINAL}")
    pred_fn = closed_neighborhood_jaccard_predicate(t_best)
    xs = [round(s / steps, 6) for s in range(steps + 1)]
    sum_scores = [0.0] * (steps + 1)

    for run_i in range(1, AVG_FINAL + 1):
        t_run0 = time.perf_counter()
        G = G_base.copy()
        edge_stack = list(G.edges())
        rng.shuffle(edge_stack)

        for s in range(steps):
            pred = get_communities(G, pred_fn)
            pred_pair_count = build_pair_counts(pred)
            score = omega_index(
                ground_truth_pair_counts=gt_pair_counts,
                detected_pair_counts=pred_pair_count,
                total_pairs=total_pairs,
            )
            sum_scores[s] += score
            rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

        pred = get_communities(G, pred_fn)
        pred_pair_count = build_pair_counts(pred)
        score = omega_index(
            ground_truth_pair_counts=gt_pair_counts,
            detected_pair_counts=pred_pair_count,
            total_pairs=total_pairs,
        )
        sum_scores[steps] += score

        dt_run = time.perf_counter() - t_run0
        print(f"      [RUN {run_i}/{AVG_FINAL}] dt={dt_run:.2f}s")

    ys = [v / AVG_FINAL for v in sum_scores]

    for p, score in zip(xs, ys):
        write_sweep_row(sweep_writer, overlap, n, z, ring_size, t_best, p, score)

    p_crit = None
    for x, y in zip(xs, ys):
        if y < P_CRIT_THRESH:
            p_crit = x
            break

    write_summary_row(summary_writer, overlap, n, z, ring_size, t_best, p_crit)
    print(f"    [RESULT] t_best={t_best:.2f}  p_crit={p_crit}")

    return t_best, p_crit


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_experiment(
    out_dir: str = "results/pcrit_experiment",
    seed: int = 42,
):
    configs = remaining_configs()

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    root_dir = os.path.join(out_dir, f"run_{ts}")
    os.makedirs(root_dir, exist_ok=True)

    search_path  = os.path.join(root_dir, "threshold_search.csv")
    sweep_path   = os.path.join(root_dir, "p_sweep.csv")
    summary_path = os.path.join(root_dir, "summary.csv")

    search_fields  = ["overlap", "n", "z", "ring_size", "phase", "threshold", "p_crit_est"]
    sweep_fields   = ["overlap", "n", "z", "ring_size", "t_best", "p", "score"]
    summary_fields = ["overlap", "n", "z", "ring_size", "t_best", "p_crit"]

    rng = random.Random(seed)

    total = len(configs)
    t_global0 = time.perf_counter()

    print(f"Total configs: {total}")
    print(f"Output: {root_dir}\n")

    with (
        open(search_path,  "w", newline="") as sf,
        open(sweep_path,   "w", newline="") as pf,
        open(summary_path, "w", newline="") as smf,
    ):
        search_writer  = csv.DictWriter(sf,  fieldnames=search_fields);  search_writer.writeheader()
        sweep_writer   = csv.DictWriter(pf,  fieldnames=sweep_fields);   sweep_writer.writeheader()
        summary_writer = csv.DictWriter(smf, fieldnames=summary_fields); summary_writer.writeheader()

        for ci, (overlap, n, z, ring_size) in enumerate(configs, start=1):
            t_cfg0 = time.perf_counter()
            print(f"\n[CONFIG {ci}/{total}] overlap={overlap}, n={n}, z={z}, ring_size={ring_size}")

            run_config(
                overlap=overlap, n=n, z=z, ring_size=ring_size,
                rng=rng,
                search_writer=search_writer,
                sweep_writer=sweep_writer,
                summary_writer=summary_writer,
            )

            sf.flush()
            pf.flush()
            smf.flush()

            dt_cfg = time.perf_counter() - t_cfg0
            elapsed = time.perf_counter() - t_global0
            remaining = (elapsed / ci) * (total - ci)
            print(f"[CONFIG {ci}/{total}] done  dt={dt_cfg:.1f}s  elapsed={elapsed:.1f}s  est_remaining={remaining:.1f}s")

    dt_all = time.perf_counter() - t_global0
    print(f"\n[ALL DONE] total dt={dt_all:.1f}s")
    print(f"Results in: {root_dir}")


if __name__ == "__main__":
    run_experiment()
