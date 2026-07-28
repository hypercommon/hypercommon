"""
generate_dataset — produce the GNN training corpus for hypercommon t_best prediction.

For each (shape, overlap_pct, run) trajectory:
  1. Build ring_lattice(sizes, zs) and apply random overlap merge.
  2. Walk p = 0..1 in steps of P_STEP. At each p, sweep t over T_GRID and record
     omega(p, t). Then rewire k_step edges and continue.
  3. Save the initial graph + omega grid + per-p argmax labels.

Each trajectory writes a 'done' marker file as its very last step. Re-running
the script SKIPS any trajectory that already has 'done', so it picks up exactly
where it left off across restarts / power loss / kills. Partial output from a
crashed trajectory is overwritten on its next run.

Output layout:
  results/gnn/rewire_random/dataset_v2/
    manifest.csv                  rebuilt on start by scanning 'done' markers
    progress.log                  append-only timing/status log
    trajectories/<traj_id>/
      init.edges.gz               edge list of G after overlap merge, BEFORE rewiring (==p_000)
      init_meta.json              {shape, sizes, zs, n, rings, overlap_pct, run, seed,
                                   M_actual, k_step, n_actual, merged_pairs}
      graphs.npz                  101 keys 'p_000'..'p_100', each (E, 2) int32 edge array
      omega_grid.npz              (n_p, n_t) float array, omega[p_idx, t_idx]
      labels.csv                  per-p: p, t_argmax, omega_at_argmax
      done                        written last

Trajectory id format:  {shape_name}_o{overlap_pct:02d}_run{run}
  (overlap_pct is the integer percent, so o00..o10)

dataset_v1 holds the earlier uniform-ring corpus, whose ids and meta use the
(n, z, ring) schema this script no longer emits; it is left untouched.

Run from venv:
  ./.venv/Scripts/python.exe -m experiments.gnn.generate_dataset
  HC_ONLY_SHAPES=u4x100z8 ./.venv/Scripts/python.exe -m experiments.gnn.generate_dataset
"""

from __future__ import annotations

import os
import csv
import json
import gzip
import glob
import time
import random
import hashlib
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from fractions import Fraction

import numpy as np

warnings.filterwarnings("ignore")

from generators.overlap import apply_overlap, ground_truth_with_overlap
from generators.ring_lattice import ring_lattice, ring_lattice_edge_count
from utils.rewiring import rewire_step
from metrics.omega import build_pair_counts


# =====================================================================
# Config — what to sweep
# =====================================================================

# A shape is (name, sizes, zs): ring r has sizes[r] nodes at degree zs[r], so
# n = sum(sizes) and the ring count is len(sizes). Listing shapes explicitly
# rather than taking a cartesian product of (n, z, ring_size) is what allows
# rings to differ from one another within a single graph.
#
# Names appear in trajectory ids and output paths, so keep them short, unique,
# and stable — renaming one orphans everything already computed under it.
SHAPES: list[tuple[str, list[int], list[int]]] = [
    # --- uniform baselines, equivalent to the old (n, z, ring_size) grid ---
    ("u4x100z8",   [100] * 4, [8] * 4),
    ("u4x100z16",  [100] * 4, [16] * 4),
]

# overlap as integer percent so traj IDs are clean ints  (0..10 -> 0.00..0.10)
OVERLAP_PCTS = list(range(0, 11))     # 0, 1, ..., 10  (i.e., 0% .. 10%)

RUNS_PER_CONFIG = 5

P_STEP = 0.01     # p in {0, 0.01, ..., 1.00}  -> 101 values
T_STEP = 0.01     # t in {0, 0.01, ..., 1.00}  -> 101 values

# Parallel workers for the t-sweep. Override with HC_WORKERS to cap CPU usage
# (the box this was developed on has 22 logical cores; 7 is ~32%).
N_WORKERS = int(os.environ.get("HC_WORKERS", "7"))

# Optional: restrict this run to named shapes, e.g. HC_ONLY_SHAPES=u4x100z8,u4x100z16
# to finish one group before starting the next. Empty means "all of SHAPES".
_ONLY_SHAPES = {x.strip() for x in os.environ.get("HC_ONLY_SHAPES", "").split(",") if x.strip()}

OUTPUT_ROOT = os.path.join("results", "gnn", "rewire_random", "dataset_v2")


# =====================================================================
# Helpers
# =====================================================================

def t_grid() -> list[float]:
    n_pts = int(round(1.0 / T_STEP)) + 1
    return [round(i * T_STEP, 4) for i in range(n_pts)]


def trajectory_id(shape_name: str, overlap_pct: int, run: int) -> str:
    return f"{shape_name}_o{overlap_pct:02d}_run{run}"


def trajectory_seed(traj_id: str) -> int:
    """Deterministic seed derived from the trajectory id, so reruns are idempotent."""
    h = hashlib.md5(traj_id.encode()).hexdigest()
    return int(h[:8], 16)


def all_trajectories() -> list[dict]:
    """All (shape, overlap_pct, run) combinations in canonical iteration order."""
    out = []
    for name, sizes, zs in SHAPES:
        for op in OVERLAP_PCTS:
            for run in range(1, RUNS_PER_CONFIG + 1):
                tid = trajectory_id(name, op, run)
                out.append({
                    "id":          tid,
                    "shape":       name,
                    "sizes":       list(sizes),
                    "zs":          list(zs),
                    "n":           sum(sizes),
                    "rings":       len(sizes),
                    "overlap_pct": op,
                    "overlap":     op / 100.0,
                    "run":         run,
                    "seed":        trajectory_seed(tid),
                })
    return out


def validate_rewiring_plan(M0: int, p_step: float) -> int:
    p = Fraction(str(p_step)).limit_denominator()
    if p.numerator != 1:
        raise ValueError(f"p_step must be 1/steps. Got {p_step}")
    steps = p.denominator
    if M0 % steps != 0:
        raise ValueError(f"M0={M0} not divisible by steps={steps}")
    return steps


# =====================================================================
# Worker — one (graph snapshot, t) -> omega
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
# One trajectory
# =====================================================================

def write_init_graph_gz(path: str, edges: list[tuple[int, int]]) -> None:
    with gzip.open(path, "wt", newline="") as f:
        for u, v in edges:
            f.write(f"{u} {v}\n")


def run_trajectory(traj: dict, root: str, pool: ProcessPoolExecutor) -> dict:
    """Compute one trajectory and write all its files. Returns manifest row dict."""
    traj_dir = os.path.join(root, "trajectories", traj["id"])
    os.makedirs(traj_dir, exist_ok=True)

    # Wipe any partial files from a previous crashed attempt
    for fname in ("init.edges.gz", "init_meta.json", "graphs.npz", "omega_grid.npz", "labels.csv", "done"):
        p = os.path.join(traj_dir, fname)
        if os.path.exists(p):
            os.remove(p)

    rng = random.Random(traj["seed"])

    # Build initial graph and ground truth
    sizes, zs = traj["sizes"], traj["zs"]
    G = ring_lattice(sizes, zs)
    merged = apply_overlap(G, sizes, overlap=traj["overlap"], rng=rng)
    truth = ground_truth_with_overlap(sizes, merged)
    gt_pc = build_pair_counts(truth)
    n_actual = G.number_of_nodes()
    total_pairs = n_actual * (n_actual - 1) // 2
    M_actual = G.number_of_edges()

    M0 = ring_lattice_edge_count(sizes, zs)
    steps = validate_rewiring_plan(M0, P_STEP)
    k_step = M_actual // steps

    # Save init graph + meta
    init_edges = list(G.edges())
    write_init_graph_gz(os.path.join(traj_dir, "init.edges.gz"), init_edges)
    meta = {
        "id":           traj["id"],
        "shape":        traj["shape"],
        "sizes":        sizes,
        "zs":           zs,
        "n":            traj["n"],
        "rings":        traj["rings"],
        "overlap_pct":  traj["overlap_pct"],
        "overlap":      traj["overlap"],
        "run":          traj["run"],
        "seed":         traj["seed"],
        "n_actual":     n_actual,
        "M_actual":     M_actual,
        "k_step":       k_step,
        "p_step":       P_STEP,
        "t_step":       T_STEP,
        "n_p":          steps + 1,
        "n_t":          len(t_grid()),
        "merged_pairs": [[v, u] for v, u in merged.items()],
    }
    with open(os.path.join(traj_dir, "init_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    # Walk p; at each p sweep t in parallel
    edge_stack = list(G.edges())
    rng.shuffle(edge_stack)

    TGRID = t_grid()
    n_p = steps + 1
    n_t = len(TGRID)
    omega_grid = np.full((n_p, n_t), np.nan, dtype=np.float64)
    graphs_dict: dict[str, np.ndarray] = {}

    labels_rows = []

    t_traj0 = time.perf_counter()
    for s in range(n_p):
        edges = list(G.edges())
        nodes = list(G.nodes())

        # Snapshot the graph at this p step
        graphs_dict[f"p_{s:03d}"] = np.array(edges, dtype=np.int32)

        futures = {t: pool.submit(_hypercommon_omega_at_t, t, edges, nodes, gt_pc, total_pairs)
                   for t in TGRID}

        row = np.full(n_t, np.nan, dtype=np.float64)
        for ti, t in enumerate(TGRID):
            _, omega = futures[t].result()
            row[ti] = omega
        omega_grid[s, :] = row

        # argmax + omega_at_argmax
        valid_mask = ~np.isnan(row)
        if valid_mask.any():
            ti_arg = int(np.argmax(np.where(valid_mask, row, -np.inf)))
            t_arg = TGRID[ti_arg]
            om_arg = float(row[ti_arg])
        else:
            t_arg = float("nan")
            om_arg = float("nan")
        labels_rows.append({"p": round(s / steps, 6), "t_argmax": t_arg, "omega_at_argmax": om_arg})

        if s < steps:
            rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

    # Persist
    np.savez_compressed(os.path.join(traj_dir, "graphs.npz"), **graphs_dict)
    np.savez_compressed(os.path.join(traj_dir, "omega_grid.npz"),
                        omega=omega_grid,
                        p_grid=np.array([round(i / steps, 6) for i in range(n_p)]),
                        t_grid=np.array(TGRID))

    with open(os.path.join(traj_dir, "labels.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["p", "t_argmax", "omega_at_argmax"])
        w.writeheader()
        for row in labels_rows:
            w.writerow(row)

    # Mark complete LAST
    with open(os.path.join(traj_dir, "done"), "w") as f:
        f.write(datetime.now().isoformat() + "\n")

    dt = time.perf_counter() - t_traj0

    return {
        "id":            traj["id"],
        "shape":         traj["shape"],
        "sizes":         " ".join(map(str, sizes)),
        "zs":            " ".join(map(str, zs)),
        "n":             traj["n"],
        "rings":         traj["rings"],
        "overlap_pct":   traj["overlap_pct"],
        "run":           traj["run"],
        "seed":          traj["seed"],
        "n_actual":      n_actual,
        "M_actual":      M_actual,
        "k_step":        k_step,
        "n_p":           n_p,
        "n_t":           n_t,
        "duration_sec":  round(dt, 2),
        "completed_at":  datetime.now().isoformat(),
    }


# =====================================================================
# Resume: scan completed
# =====================================================================

def scan_completed(root: str) -> set[str]:
    out = set()
    for d in glob.glob(os.path.join(root, "trajectories", "*")):
        if os.path.isdir(d) and os.path.exists(os.path.join(d, "done")):
            out.add(os.path.basename(d))
    return out


def rebuild_manifest_from_disk(root: str) -> list[dict]:
    rows = []
    for d in sorted(glob.glob(os.path.join(root, "trajectories", "*"))):
        if not os.path.exists(os.path.join(d, "done")):
            continue
        meta_path = os.path.join(d, "init_meta.json")
        if not os.path.exists(meta_path):
            continue
        meta = json.load(open(meta_path))
        sizes = meta.get("sizes") or []
        zs = meta.get("zs") or []
        rows.append({
            "id":          meta.get("id"),
            "shape":       meta.get("shape"),
            "sizes":       " ".join(map(str, sizes)),
            "zs":          " ".join(map(str, zs)),
            "n":           meta.get("n"),
            "rings":       meta.get("rings"),
            "overlap_pct": meta.get("overlap_pct"),
            "run":         meta.get("run"),
            "seed":        meta.get("seed"),
            "n_actual":    meta.get("n_actual"),
            "M_actual":    meta.get("M_actual"),
            "k_step":      meta.get("k_step"),
            "n_p":         meta.get("n_p"),
            "n_t":         meta.get("n_t"),
        })
    return rows


MANIFEST_FIELDS = ["id", "shape", "sizes", "zs", "n", "rings", "overlap_pct", "run",
                   "seed", "n_actual", "M_actual", "k_step", "n_p", "n_t"]


def write_manifest(root: str, rows: list[dict]) -> None:
    path = os.path.join(root, "manifest.csv")
    fields = MANIFEST_FIELDS
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


# =====================================================================
# Main
# =====================================================================

def run(out_root: str = OUTPUT_ROOT):
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(os.path.join(out_root, "trajectories"), exist_ok=True)

    log_path = os.path.join(out_root, "progress.log")

    def log(msg: str):
        line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        print(line, flush=True)
        with open(log_path, "a") as f:
            f.write(line + "\n")

    all_trajs = all_trajectories()
    completed = scan_completed(out_root)

    # Rebuild manifest from what's actually on disk (handles missing/corrupt manifest)
    write_manifest(out_root, rebuild_manifest_from_disk(out_root))

    pending = [t for t in all_trajs if t["id"] not in completed]
    if _ONLY_SHAPES:
        unknown = _ONLY_SHAPES - {name for name, _, _ in SHAPES}
        if unknown:
            raise SystemExit(f"HC_ONLY_SHAPES names no such shape: {sorted(unknown)}")
        skipped = [t for t in pending if t["shape"] not in _ONLY_SHAPES]
        pending = [t for t in pending if t["shape"] in _ONLY_SHAPES]

    log(f"=== generate_dataset start ===")
    log(f"output: {out_root}")
    log(f"total trajectories: {len(all_trajs)}, completed: {len(completed)}, pending: {len(pending)}")
    log(f"workers: {N_WORKERS}")
    if _ONLY_SHAPES:
        log(f"restricted to shapes {sorted(_ONLY_SHAPES)} "
            f"({len(skipped)} pending trajectories deferred)")

    if not pending:
        log("nothing to do; everything already complete.")
        return

    t_global0 = time.perf_counter()

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for i, traj in enumerate(pending, start=1):
            try:
                row = run_trajectory(traj, out_root, pool)
            except Exception as e:
                log(f"FAILED {traj['id']}: {type(e).__name__}: {e}")
                continue

            elapsed = time.perf_counter() - t_global0
            avg_dt = elapsed / i
            est_remaining = avg_dt * (len(pending) - i)
            log(f"  done {traj['id']}  dt={row['duration_sec']:.1f}s  "
                f"({i}/{len(pending)})  elapsed={elapsed:.0f}s  est_rem={est_remaining:.0f}s")

            # Append to manifest after every completed trajectory
            with open(os.path.join(out_root, "manifest.csv"), "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
                w.writerow({k: row.get(k) for k in MANIFEST_FIELDS})

    log(f"=== generate_dataset done ===  total dt={time.perf_counter()-t_global0:.0f}s")


if __name__ == "__main__":
    run()
