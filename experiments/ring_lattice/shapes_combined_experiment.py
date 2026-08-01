"""
shapes_combined_experiment - hypercommon vs 4 baselines when SIZE and DEGREE both vary.

The two earlier sweeps each moved one factor:

  shapes_algos_experiment   ring size uniform, z varies    -> hypercommon won all 5 shapes
  shapes_sizes_experiment   z uniform, ring size varies    -> leiden won or tied on 5 of 8

Read together they locate the boundary: hypercommon leads on sparse, uneven
graphs and trails on dense, even ones. The p_crit gap against leiden at
overlap=0 tracked mean degree and size ratio almost monotonically:

  mean degree  size ratio   gap
      8           5.0      +0.32
      8           4.0      +0.17
      8           2.0       0.00
     16           4.0      +0.02
     16           2.0      -0.03
     32           2.0      -0.05

This sweep moves both factors at once, which neither earlier one could do, and
so can ask questions they could not:

  c_grad_aligned vs c_grad_inverse   identical sizes; z aligned with size
                                     (big rings dense) or inverted (big rings
                                     sparse). Does the size-degree CORRELATION
                                     matter, independent of either factor?
  c_smallhub vs c_bigsparse          identical partition, degrees swapped
                                     between the small and large rings.
  c_dense_mild vs c_dense_uneven     both dense; does leiden's advantage there
                                     survive size unevenness?

Every shape varies BOTH size and z - single-factor shapes belong in the two
earlier scripts.

For each (shape, overlap%, run):
  1. Build ring_lattice(sizes, zs) and merge floor(overlap * n) inter-ring node
     pairs, so the ground truth genuinely overlaps.
  2. Walk p = 0..1 in P_STEP increments, rewiring k_step edges per step.
  3. At every p, run all algorithms on the same graph snapshot, in parallel, and
     score each against the ground truth with the omega index.

Hypercommon is evaluated at its best threshold per snapshot rather than one fixed
t: a coarse grid at p=0, then an adaptive window around the previous step's
winner. Both the winning t and the searched range are recorded, so the plotting
stage can show how t* moves without re-running anything.

Output - results/shapes_combined/run_<ts>/
  records.csv     one row per (shape, overlap, run, p, algo); flushed as the
                  sweep proceeds, so a killed run leaves usable data
  shapes.csv      one row per shape: sizes, zs, n, rings, edges, mean degree
  progress.log    append-only timing and status
  done/<key>      marker per completed (shape, overlap, run); re-running skips
                  any combination that already has one

records.csv columns
  shape run_id overlap_pct run p step
  n rings ring_size zs n_actual edges_actual k_step
  algo omega n_communities elapsed_sec
  t_argmax t_grid_lo t_grid_hi t_grid_size

ring_size holds the FIRST ring's size, since sizes differ within a shape; the
full list is in the sizes column of shapes.csv. This script only computes and
writes; plotting is a separate step over records.csv.

Concurrency note: eight simultaneous shapes exhausted memory during the sizes
sweep and killed two of them mid-run. Launch at most four or five at a time.

Run from venv:
  ./.venv/Scripts/python.exe -m experiments.ring_lattice.shapes_combined_experiment
  ./.venv/Scripts/python.exe -m experiments.ring_lattice.shapes_combined_experiment --validate --build
  HC_ONLY_SHAPES=c_skew_sparse ./.venv/Scripts/python.exe -m experiments.ring_lattice.shapes_combined_experiment
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import random
import time
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from fractions import Fraction

# cdlib initialises global state on import; do it in the parent before any
# worker process forks off.
warnings.filterwarnings("ignore")
from cdlib import algorithms as _cdlib_preload  # noqa: F401

from generators.overlap import apply_overlap, ground_truth_with_overlap
from generators.ring_lattice import (
    ring_communities,
    ring_lattice,
    ring_lattice_edge_count,
)
from metrics.omega import build_pair_counts
from utils.rewiring import rewire_step


# =====================================================================
# Config
# =====================================================================

# Ordered MOST EXPENSIVE FIRST, so the long poles start immediately rather than
# becoming a tail once the cheap shapes finish.
#
# Every shape varies BOTH ring size and z. Constraints: n = 2000, every ring
# needs z < its own size, and the edge count must be divisible by 100 for p to
# reach exactly 1.0.
SHAPES: list[tuple[str, list[int], list[int]]] = [
    # --- dense: leiden's favourable regime, now with uneven sizes ---
    ("c_dense_mild",
     [300, 300, 250, 250, 200, 200, 200, 100, 100, 100],
     [40, 32, 32, 32, 24, 24, 24, 16, 16, 16]),
    ("c_dense_uneven",
     [400, 300, 200, 200, 200, 150, 150, 150, 150, 100],
     [40, 32, 32, 24, 24, 24, 16, 16, 16, 16]),

    # two disjoint regimes in one graph: 3 big dense rings, 16 small sparse ones
    ("c_two_worlds",
     [400] * 3 + [50] * 16,
     [32] * 3 + [8] * 16),

    # --- the size-degree correlation pair: same partition, degrees swapped ---
    ("c_smallhub",  [50] * 20 + [200] * 5, [8] * 20 + [32] * 5),
    ("c_bigsparse", [50] * 20 + [200] * 5, [24] * 20 + [8] * 5),

    # --- the alignment pair: identical sizes, z aligned then inverted ---
    ("c_grad_aligned",
     [400, 300, 250, 200, 200, 150, 150, 100, 100, 150],
     [32, 24, 24, 16, 16, 16, 12, 12, 8, 8]),
    ("c_grad_inverse",
     [400, 300, 250, 200, 200, 150, 150, 100, 100, 150],
     [8, 8, 12, 12, 16, 16, 24, 24, 32, 32]),

    # --- sparse and uneven: where hypercommon led by the widest margin ---
    ("c_mixed4_hetz",
     [50] * 5 + [100] * 5 + [200] * 5 + [250],
     [8] * 5 + [12] * 5 + [16] * 5 + [24]),
    ("c_skew_sparse",
     [500, 300, 200, 150, 150, 150, 150, 100, 100, 200],
     [8, 8, 8, 12, 12, 16, 16, 16, 24, 24]),
    # the most uneven partition in any sweep: a 20x size ratio
    ("c_sparse_extreme",
     [600, 400, 300, 200, 150, 120, 100, 60, 40, 30],
     [16, 12, 12, 8, 8, 8, 8, 8, 8, 8]),
]

# Overlap as integer percent, so it stays exact in ids and output rows.
OVERLAP_PCTS = [0, 2, 5, 10]

RUNS_PER_CONFIG = 5

# slpa and label_propagation are omitted: on the completed uniform runs neither
# ever exceeded omega 0.50, including at p=0 where every algorithm below scores
# 1.000. They cannot lose structure they never recovered, so they add rows
# without adding a comparison.
ALGOS = [
    "hypercommon",
    "leiden",
    "walktrap",
    "demon",
    "angel",
]

P_STEP = 0.01

# Hypercommon threshold search: a coarse sweep at p=0 to find the peak, then a
# narrow window around the previous winner, since t* moves slowly in p.
T_INIT_STEP = 0.04
T_ADAPTIVE_HALF = 0.04
T_ADAPTIVE_STEP = 0.01

N_WORKERS = int(os.environ.get("HC_WORKERS", str(len(ALGOS))))

_ONLY_SHAPES = {x.strip() for x in os.environ.get("HC_ONLY_SHAPES", "").split(",") if x.strip()}

OUTPUT_ROOT = os.path.join("results", "shapes_combined")

# Only len(ALGOS) tasks exist per p-step, so a bigger pool cannot speed up one
# shape. To use more of the machine, launch several shapes at once with
# HC_ONLY_SHAPES and point them at a shared HC_OUT_DIR: each writes its own
# records part file, so there is no contention.
_OUT_DIR = os.environ.get("HC_OUT_DIR", "").strip()

# Distinguishes concurrent writers inside one output directory.
_RUN_TAG = os.environ.get("HC_RUN_TAG", "").strip()

RECORD_FIELDS = [
    "shape", "run_id", "overlap_pct", "run", "p", "step",
    "n", "rings", "ring_size", "zs", "n_actual", "edges_actual", "k_step",
    "algo", "omega", "n_communities", "elapsed_sec",
    "t_argmax", "t_grid_lo", "t_grid_hi", "t_grid_size",
]

SHAPE_FIELDS = [
    "shape", "sizes", "zs", "n", "rings", "ring_size",
    "edges", "mean_degree", "distinct_z", "z_min", "z_max",
]


# =====================================================================
# Planning
# =====================================================================

def validate_rewiring_plan(edges: int, p_step: float) -> int:
    """Number of rewiring steps for p to reach exactly 1.0."""
    p = Fraction(str(p_step)).limit_denominator()
    if p.numerator != 1:
        raise ValueError(f"p_step must be 1/steps. Got {p_step}")
    steps = p.denominator
    if edges % steps != 0:
        raise ValueError(f"edge count {edges} not divisible by steps={steps}")
    return steps


def run_key(shape: str, overlap_pct: int, run: int) -> str:
    return f"{shape}_o{overlap_pct:02d}_run{run}"


def all_units() -> list[dict]:
    """Every (shape, overlap, run) the sweep covers, in canonical order."""
    units = []
    for name, sizes, zs in SHAPES:
        for overlap_pct in OVERLAP_PCTS:
            for run in range(1, RUNS_PER_CONFIG + 1):
                units.append({
                    "shape": name,
                    "sizes": list(sizes),
                    "zs": list(zs),
                    "overlap_pct": overlap_pct,
                    "run": run,
                    "run_id": run_key(name, overlap_pct, run),
                })
    return units


def unit_seed(run_id: str) -> int:
    """Deterministic per-unit seed so a resumed run reproduces the same graphs."""
    return int(hashlib.md5(run_id.encode()).hexdigest()[:8], 16)


def make_t_grid(step: int, prev_t_best: float | None) -> list[float]:
    """Coarse grid at p=0; a window around the previous winner afterwards."""
    if step == 0 or prev_t_best is None:
        count = int(round((1.0 - T_INIT_STEP) / T_INIT_STEP)) + 1
        return [round((i + 1) * T_INIT_STEP, 4)
                for i in range(count)
                if (i + 1) * T_INIT_STEP < 1.0]

    half = int(round(T_ADAPTIVE_HALF / T_ADAPTIVE_STEP))
    grid = []
    for offset in range(-half, half + 1):
        t = round(prev_t_best + offset * T_ADAPTIVE_STEP, 4)
        if 0.01 <= t <= 0.99:
            grid.append(t)
    return grid


# =====================================================================
# Worker — one algorithm on one graph snapshot
# =====================================================================

def _run_algo(algo, params, edges, nodes, gt_pair_counts, total_pairs):
    """Score one algorithm on one snapshot.

    Returns (algo, omega, n_communities, elapsed_sec, t_argmax_or_None). Failures
    come back as NaN rather than raising, so one bad algorithm cannot abort a
    sweep that has been running for days.
    """
    import warnings as _warnings
    _warnings.filterwarnings("ignore")

    import time as _time
    import networkx as nx

    from metrics.omega import build_pair_counts as _pair_counts, omega_index as _omega

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    started = _time.perf_counter()

    if algo == "hypercommon":
        from hypercommon.algorithm import get_communities as _detect
        from predicates.jaccard import closed_neighborhood_jaccard_predicate as _predicate

        best_omega = float("-inf")
        best_t = None
        best_count = 0
        for t in params["t_grid"]:
            try:
                communities = _detect(G, _predicate(t))
                score = _omega(gt_pair_counts, _pair_counts(communities), total_pairs)
            except Exception:
                continue
            if score > best_omega:
                best_omega, best_t, best_count = score, t, len(communities)

        elapsed = _time.perf_counter() - started
        if best_t is None:
            return algo, float("nan"), 0, elapsed, None
        return algo, float(best_omega), best_count, elapsed, float(best_t)

    try:
        from cdlib import algorithms as A

        if algo == "leiden":
            found = A.leiden(G)
        elif algo == "label_propagation":
            found = A.label_propagation(G)
        elif algo == "walktrap":
            found = A.walktrap(G)
        elif algo == "slpa":
            found = A.slpa(G)
        elif algo == "demon":
            found = A.demon(G, epsilon=0.25)
        elif algo == "angel":
            found = A.angel(G, threshold=0.25)
        else:
            raise ValueError(f"unknown algo: {algo}")

        communities = [set(c) for c in found.communities]
        score = _omega(gt_pair_counts, _pair_counts(communities), total_pairs)
        return algo, float(score), len(communities), _time.perf_counter() - started, None
    except Exception:
        return algo, float("nan"), 0, _time.perf_counter() - started, None


# =====================================================================
# One unit — one (shape, overlap, run) walked across p
# =====================================================================

def run_unit(unit: dict, writer: csv.DictWriter, handle, pool: ProcessPoolExecutor, log) -> None:
    sizes, zs = unit["sizes"], unit["zs"]
    overlap = unit["overlap_pct"] / 100.0

    rng = random.Random(unit_seed(unit["run_id"]))

    G = ring_lattice(sizes, zs)
    merged = apply_overlap(G, sizes, overlap, rng)
    truth = ground_truth_with_overlap(sizes, merged)
    gt_pair_counts = build_pair_counts(truth)

    n_actual = G.number_of_nodes()
    total_pairs = n_actual * (n_actual - 1) // 2
    edges_actual = G.number_of_edges()

    steps = validate_rewiring_plan(ring_lattice_edge_count(sizes, zs), P_STEP)
    k_step = edges_actual // steps

    edge_stack = list(G.edges())
    rng.shuffle(edge_stack)

    base = {
        "shape": unit["shape"],
        "run_id": unit["run_id"],
        "overlap_pct": unit["overlap_pct"],
        "run": unit["run"],
        "n": sum(sizes),
        "rings": len(sizes),
        "ring_size": sizes[0],
        "zs": " ".join(map(str, zs)),
        "n_actual": n_actual,
        "edges_actual": edges_actual,
        "k_step": k_step,
    }

    prev_t_best = None
    started = time.perf_counter()

    for step in range(steps + 1):
        snapshot_edges = list(G.edges())
        snapshot_nodes = list(G.nodes())
        t_grid = make_t_grid(step, prev_t_best)

        futures = {
            algo: pool.submit(
                _run_algo,
                algo,
                {"t_grid": t_grid} if algo == "hypercommon" else {},
                snapshot_edges, snapshot_nodes, gt_pair_counts, total_pairs,
            )
            for algo in ALGOS
        }

        for algo in ALGOS:
            name, omega, n_communities, elapsed, t_argmax = futures[algo].result()

            row = dict(base)
            row.update({
                "p": round(step / steps, 6),
                "step": step,
                "algo": name,
                "omega": omega,
                "n_communities": n_communities,
                "elapsed_sec": round(elapsed, 4),
                "t_argmax": t_argmax if t_argmax is not None else "",
                "t_grid_lo": min(t_grid) if name == "hypercommon" else "",
                "t_grid_hi": max(t_grid) if name == "hypercommon" else "",
                "t_grid_size": len(t_grid) if name == "hypercommon" else "",
            })
            writer.writerow(row)

            if name == "hypercommon" and t_argmax is not None:
                prev_t_best = t_argmax

        handle.flush()

        if step < steps:
            rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

    log(f"  done {unit['run_id']}  n={sum(sizes)} edges={edges_actual} "
        f"merged={len(merged)} dt={time.perf_counter() - started:.1f}s")


def write_shapes_table(path: str) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SHAPE_FIELDS)
        writer.writeheader()
        for name, sizes, zs in SHAPES:
            edges = ring_lattice_edge_count(sizes, zs)
            writer.writerow({
                "shape": name,
                "sizes": " ".join(map(str, sizes)),
                "zs": " ".join(map(str, zs)),
                "n": sum(sizes),
                "rings": len(sizes),
                "ring_size": sizes[0],
                "edges": edges,
                "mean_degree": round(2 * edges / sum(sizes), 3),
                "distinct_z": len(set(zs)),
                "z_min": min(zs),
                "z_max": max(zs),
            })


# =====================================================================
# Validation — check the plan before committing days of compute
# =====================================================================

# ms per get_communities call ~ exp(a) * n^b * (E/n)^c, fitted to the measured
# calls in COST_SAMPLES. The third term is E/n, i.e. HALF the mean degree —
# feeding it the mean degree instead inflates every estimate by 2**2.35 ~ 5.1x.
COST_FIT = (-7.18, 1.24, 2.35)

COST_SAMPLES = [
    (400, 1600, 36.2), (400, 3200, 158.5), (400, 4000, 323.3),
    (800, 6400, 352.3), (1200, 9600, 600.7),
    (2000, 18400, 2120.6), (2000, 16000, 1187.8),
]

# A full unit of n=400 ring=100 z=[8,8,16,16] measured 153.6 s against
# 924 x 86.5 ms = 79.9 s of predicted bare calls. The remainder is per-call graph
# rebuilds in the worker, process dispatch and the rewiring itself.
MEASURED_OVERHEAD = 1.92


def predicted_ms(n: int, edges: int) -> float:
    a, b, c = COST_FIT
    return math.exp(a + b * math.log(n) + c * math.log(edges / n))


def measured_ms(sizes: list[int], zs: list[int]) -> float:
    """Time real hypercommon calls on this shape, near the thresholds that matter."""
    from hypercommon.algorithm import get_communities
    from predicates.jaccard import closed_neighborhood_jaccard_predicate as predicate

    G = ring_lattice(sizes, zs)
    started = time.perf_counter()
    for t in (0.10, 0.15, 0.25):
        get_communities(G, predicate(t))
    return (time.perf_counter() - started) / 3 * 1000


def check_shape(sizes: list[int], zs: list[int]) -> list[str]:
    problems = []

    if len(sizes) != len(zs):
        problems.append(f"sizes/zs length mismatch: {len(sizes)} vs {len(zs)}")
        return problems

    if len(set(sizes)) < 2:
        problems.append("ring sizes do not vary - this sweep varies BOTH factors, "
                        "so a uniform-size shape belongs in shapes_algos_experiment")

    if len(set(zs)) < 2:
        problems.append("z does not vary across rings - this sweep varies BOTH "
                        "factors, so a uniform-z shape belongs in shapes_sizes_experiment")

    for index, (size, z) in enumerate(zip(sizes, zs)):
        if z % 2:
            problems.append(f"ring {index}: z={z} is odd")
        if z >= size:
            problems.append(f"ring {index}: z={z} >= size={size}")

    try:
        validate_rewiring_plan(ring_lattice_edge_count(sizes, zs), P_STEP)
    except ValueError as exc:
        problems.append(f"rewiring plan: {exc}")

    return problems


def validate(build: bool, measure: bool) -> int:
    import networkx as nx

    steps = int(round(1 / P_STEP))
    per_shape = len(OVERLAP_PCTS) * RUNS_PER_CONFIG
    # Coarse grid at p=0, then the adaptive window at every later step.
    hc_calls = len(make_t_grid(0, None)) + steps * len(make_t_grid(1, 0.15))

    worst = max(abs(predicted_ms(n, e) / ms - 1) for n, e, ms in COST_SAMPLES)
    print(f"cost model vs {len(COST_SAMPLES)} measured calls: worst error {worst:.0%}")
    if worst > 0.5:
        print("!! cost model disagrees with its own measurements — check the units")
    print()

    print(f"algos={len(ALGOS)}: {ALGOS}")
    print(f"shapes={len(SHAPES)}  overlaps={OVERLAP_PCTS}  runs={RUNS_PER_CONFIG}")
    print(f"p steps per unit: {steps + 1}   units per shape: {per_shape}")
    print(f"hypercommon threshold evaluations per unit: {hc_calls}"
          f"  (coarse {len(make_t_grid(0, None))} at p=0, "
          f"then {len(make_t_grid(1, 0.15))} per step)")
    print("runtime is hypercommon's serial t-sweep; the other algos overlap with it"
          f"{' — ms/call MEASURED per shape' if measure else ' — ms/call modelled'}")
    print()

    header = f"{'name':24} {'rings':>9} {'n':>5} {'E':>6} {'mdeg':>5} {'ms/call':>8} {'hours':>7}"
    print(header)
    print("-" * len(header))

    failures = 0
    total_hours = 0.0

    for name, sizes, zs in SHAPES:
        problems = check_shape(sizes, zs)
        n = sum(sizes)
        edges = ring_lattice_edge_count(sizes, zs)

        ms = measured_ms(sizes, zs) if measure else predicted_ms(n, edges)

        # The seven algorithms run concurrently, but hypercommon walks its t-grid
        # SEQUENTIALLY inside one worker, so a p-step costs as long as that sweep
        # and the pool gives no speedup on it. One unit is therefore about
        # hc_calls x ms of wall clock, undivided.
        hours = hc_calls * ms * MEASURED_OVERHEAD / 1000 * per_shape / 3600
        total_hours += hours

        print(f"{name:24} {f'{sizes[0]}x{len(sizes)}':>9} {n:5} {edges:6} "
              f"{2 * edges / n:5.1f} {ms:8.1f} {hours:7.1f}")

        for problem in problems:
            failures += 1
            print(f"    !! {problem}")

        if build:
            G = ring_lattice(sizes, zs)
            assert G.number_of_nodes() == n, name
            assert G.number_of_edges() == edges, name
            assert nx.number_connected_components(G) == len(sizes), name
            for community, z in zip(ring_communities(sizes), zs):
                assert {G.degree(v) for v in community} == {z}, name

    units = all_units()
    print()
    print(f"  TOTAL : {len(SHAPES):2} shapes  {len(units):4} units  "
          f"{total_hours:7.1f} h = {total_hours / 24:.1f} days")
    print(f"  records.csv will hold {len(units) * (steps + 1) * len(ALGOS):,} rows")
    print()

    names = Counter(name for name, _, _ in SHAPES)
    duplicates = [name for name, count in names.items() if count > 1]
    if duplicates:
        failures += len(duplicates)
        print(f"!! duplicate shape names: {duplicates}")

    ids = Counter(u["run_id"] for u in units)
    clashes = [i for i, count in ids.items() if count > 1]
    if clashes:
        failures += len(clashes)
        print(f"!! duplicate run ids: {clashes[:5]}")
    else:
        print(f"run ids: {len(units)} unique, e.g. {units[0]['run_id']}")

    if build:
        print("build check: every shape constructed, degrees and components verified")

    print()
    print("VALIDATION FAILED" if failures else "VALIDATION OK")
    return 1 if failures else 0


# =====================================================================
# Main
# =====================================================================

def run(out_root: str | None = None) -> None:
    if out_root is None:
        if _OUT_DIR:
            out_root = _OUT_DIR
        else:
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            out_root = os.path.join(OUTPUT_ROOT, f"run_{stamp}")

    os.makedirs(out_root, exist_ok=True)
    done_dir = os.path.join(out_root, "done")
    os.makedirs(done_dir, exist_ok=True)

    # Concurrent launches share a directory but never a file.
    suffix = f"_{_RUN_TAG}" if _RUN_TAG else ""
    log_path = os.path.join(out_root, f"progress{suffix}.log")

    def log(message: str) -> None:
        line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}"
        print(line, flush=True)
        with open(log_path, "a") as handle:
            handle.write(line + "\n")

    write_shapes_table(os.path.join(out_root, "shapes.csv"))

    units = all_units()
    if _ONLY_SHAPES:
        unknown = _ONLY_SHAPES - {name for name, _, _ in SHAPES}
        if unknown:
            raise SystemExit(f"HC_ONLY_SHAPES names no such shape: {sorted(unknown)}")
        units = [u for u in units if u["shape"] in _ONLY_SHAPES]

    pending = [u for u in units if not os.path.exists(os.path.join(done_dir, u["run_id"]))]

    log("=== shapes_combined_experiment start ===")
    log(f"output: {out_root}")
    log(f"algos: {ALGOS}")
    log(f"shapes: {len(SHAPES)}  overlaps: {OVERLAP_PCTS}  runs: {RUNS_PER_CONFIG}")
    log(f"units: {len(units)} selected, {len(pending)} pending")
    log(f"workers: {N_WORKERS}")
    if _ONLY_SHAPES:
        log(f"restricted to shapes {sorted(_ONLY_SHAPES)}")

    if not pending:
        log("nothing to do; everything already complete.")
        return

    records_path = os.path.join(out_root, f"records{suffix}.csv")
    fresh = not os.path.exists(records_path)

    started = time.perf_counter()

    with open(records_path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_FIELDS)
        if fresh:
            writer.writeheader()
            handle.flush()

        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            for index, unit in enumerate(pending, start=1):
                try:
                    run_unit(unit, writer, handle, pool, log)
                except Exception as exc:
                    log(f"  FAILED {unit['run_id']}: {type(exc).__name__}: {exc}")
                    continue

                # Marker written last, so an interrupted unit is redone in full.
                with open(os.path.join(done_dir, unit["run_id"]), "w") as marker:
                    marker.write(datetime.now().isoformat() + "\n")

                elapsed = time.perf_counter() - started
                remaining = elapsed / index * (len(pending) - index)
                log(f"  [{index}/{len(pending)}] elapsed={elapsed:.0f}s est_rem={remaining:.0f}s")

    log(f"=== shapes_combined_experiment done ===  total dt={time.perf_counter() - started:.0f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true",
                        help="print the plan and check every shape; run nothing")
    parser.add_argument("--build", action="store_true",
                        help="with --validate: construct each graph and assert its structure")
    parser.add_argument("--measure", action="store_true",
                        help="with --validate: time real calls instead of using the cost model")
    args = parser.parse_args()

    if args.validate:
        raise SystemExit(validate(build=args.build, measure=args.measure))

    run()


if __name__ == "__main__":
    main()
