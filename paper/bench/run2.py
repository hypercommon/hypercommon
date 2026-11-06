"""
Benchmark main runs: the method set and grids of methods2.py,
and every cover scored with

  omega, onmi, f1, coverage, n_comms   as in run.py
  eq          overlapping modularity of the cover (Shen et al. 2009), the label-free
              criterion: each method's label-free choice is its parameter with the
              largest eq on that graph (analysis2.py)
  rec_mean    mean over true communities of the best Jaccard with a detected one
  rec_small / rec_large   the same over true communities no larger / larger than the
              median true size (rec_large empty when all sizes are equal)

    python run2.py <experiment> [--workers N] [--methods a,b,...] [--dry]

Rows go to results2/<experiment>.csv, appended as tasks finish; a rerun skips
tasks already present. Seeds 0..9 are test graphs, 1000..1002 calibration graphs.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from experiments import EXPERIMENTS  # noqa: E402

RESULTS = os.path.join(HERE, "results2")
FIELDS = ["experiment", "family", "p", "seed", "split", "method", "param",
          "omega", "onmi", "f1", "coverage", "n_comms", "eq", "rec_mean", "rec_small", "rec_large",
          "n_nodes", "n_edges", "n_wedges", "seconds", "error"]
TEST_SEEDS = list(range(10))
CAL_SEEDS = [1000, 1001, 1002]


def _safe(f, *a):
    """A metric that is undefined for a degenerate ground truth (one community: nulls) is NaN."""
    try:
        return f(*a)
    except ZeroDivisionError:
        return float("nan")


def run_task(args):
    exp_name, fam_id, p, seed, method = args
    import numpy as np
    import cover_metrics as M
    from hc_fast import TripleScores
    from methods2 import HC_LIKE, run_method2
    from self_tune import _adj, overlapping_modularity

    exp = EXPERIMENTS[exp_name]
    fam = exp["families"][fam_id]
    split = "cal" if seed >= 1000 else "test"
    rows = []
    t0 = time.time()
    try:
        G, truth = exp["build"](fam, p, seed)
        ts = TripleScores(G) if method in HC_LIKE else None
        A, idx = _adj(G)
        sizes = np.array([len(c) for c in truth])
        med = np.median(sizes)
        small, large = sizes <= med, sizes > med
        t_build = time.time() - t0
        t1 = time.time()
        for param, cover in run_method2(method, G, ts, seed):
            rec = M.truth_recovery(cover, truth, G.nodes())
            rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split,
                             method=method, param=param,
                             omega=M.omega(cover, truth, G.nodes()),
                             onmi=_safe(M.onmi, cover, truth, G),
                             f1=M.f1_best_match(cover, truth, G.nodes()),
                             coverage=M.coverage(cover, truth), n_comms=len(cover),
                             eq=overlapping_modularity(A, idx, cover),
                             rec_mean=float(rec.mean()),
                             rec_small=float(rec[small].mean()) if small.any() else "",
                             rec_large=float(rec[large].mean()) if large.any() else "",
                             n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(),
                             n_wedges=ts.n_wedges if ts else "",
                             seconds=round(time.time() - t1 + t_build, 3), error=""))
            t1, t_build = time.time(), 0.0
    except Exception as e:  # recorded, never swallowed silently
        rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split,
                         method=method, seconds=round(time.time() - t0, 3),
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]))
    return rows


def main():
    from methods2 import METHODS2
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--methods", default="")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    exp = EXPERIMENTS[a.experiment]
    methods = a.methods.split(",") if a.methods else (
        exp["methods"] if a.experiment == "rings_isolated" else METHODS2)
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"{a.experiment}.csv")

    done = set()
    if os.path.exists(path):
        # a task that ran out of memory is not done: it failed for want of a free machine
        # (including a module import that the memory shortage left half-initialised)
        transient = ("MemoryError", "AttributeError: module 'matplotlib'")
        with open(path, newline="") as fh:
            done = {(r["family"], f"{float(r['p']):.4f}", str(r["seed"]), r["method"])
                    for r in csv.DictReader(fh) if not (r.get("error") or "").startswith(transient)}
    tasks = []
    for fam_id in exp["families"]:
        for p in exp.get("p_grid_by_family", {}).get(fam_id, exp.get("p_grid_v2", exp["p_grid"])):
            for seed in TEST_SEEDS + CAL_SEEDS:
                for m in methods:
                    if (fam_id, f"{float(p):.4f}", str(seed), m) not in done:
                        tasks.append((a.experiment, fam_id, p, seed, m))
    slow = ["slpa", "demon", "angel", "walktrap", "leiden", "leiden_cpm"]
    tasks.sort(key=lambda t: slow.index(t[4]) if t[4] in slow else 99)
    print(f"{len(tasks)} tasks to run ({len(done)} already done) -> {path}", flush=True)
    if a.dry:
        return

    new_file = not os.path.exists(path)
    with open(path, "a", newline="") as fh, Pool(a.workers, maxtasksperchild=50) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new_file:
            w.writeheader()
        t0 = time.time()
        for i, rows in enumerate(pool.imap_unordered(run_task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            if i % 50 == 0 or i == len(tasks):
                el = time.time() - t0
                print(f"  {i}/{len(tasks)}  {el/60:.1f} min  eta {el/i*(len(tasks)-i)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
