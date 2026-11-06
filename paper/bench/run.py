"""
Benchmark runner.

    python run.py <experiment> [--workers N] [--methods a,b,...] [--dry]

One task = (family, p, seed, method). A task builds its graph deterministically,
runs the method's whole parameter grid, and scores every parameter with Omega,
ONMI and coverage. Rows are appended to results/<experiment>.csv as tasks finish,
and a rerun skips tasks already present, so an interrupted sweep resumes.

Seeds 0..9 are TEST graphs; seeds 1000..1002 are CALIBRATION graphs, used only to
freeze each method's parameter (analysis.py). The per-graph oracle is read off
the test rows.
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

RESULTS = os.path.join(HERE, "results")
FIELDS = ["experiment", "family", "p", "seed", "split", "method", "param",
          "omega", "onmi", "f1", "coverage", "n_comms", "n_nodes", "n_edges",
          "n_wedges", "seconds", "error"]
TEST_SEEDS = list(range(10))
CAL_SEEDS = [1000, 1001, 1002]


def task_key(r):
    return (r["family"], f"{float(r['p']):.4f}", str(r["seed"]), r["method"])


def run_task(args):
    exp_name, fam_id, p, seed, method = args
    import graphs
    import cover_metrics as M
    from hc_fast import TripleScores
    from methods import run_method

    exp = EXPERIMENTS[exp_name]
    fam = exp["families"][fam_id]
    split = "cal" if seed >= 1000 else "test"
    rows = []
    t0 = time.time()
    try:
        G, truth = exp["build"](fam, p, seed)
        if method.startswith("hcf_"):
            from hc_family import FamilyScores
            ts = FamilyScores(G)
        else:
            ts = TripleScores(G) if method.startswith("hc") else None
        t_build = time.time() - t0
        for param, cover in run_method(method, G, ts, seed):
            t1 = time.time()
            om = M.omega(cover, truth, G.nodes())
            on = M.onmi(cover, truth, G) if exp.get("onmi", True) else float("nan")
            rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split,
                             method=method, param=param, omega=om, onmi=on,
                             f1=M.f1_best_match(cover, truth, G.nodes()),
                             coverage=M.coverage(cover, truth), n_comms=len(cover),
                             n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(),
                             n_wedges=ts.n_wedges if ts else "", seconds=round(time.time() - t1 + t_build, 3),
                             error=""))
    except Exception as e:  # recorded, never swallowed silently
        rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split,
                         method=method, param="", omega="", onmi="", coverage="", n_comms="",
                         n_nodes="", n_edges="", n_wedges="", seconds=round(time.time() - t0, 3),
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--methods", default="")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    exp = EXPERIMENTS[a.experiment]
    methods = a.methods.split(",") if a.methods else exp["methods"]
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"{a.experiment}.csv")

    done = set()
    if os.path.exists(path):
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                done.add(task_key(r))

    tasks = []
    for fam_id in exp["families"]:
        for p in exp.get("p_grid_by_family", {}).get(fam_id, exp["p_grid"]):
            for seed in TEST_SEEDS[: exp.get("n_test", 10)] + CAL_SEEDS[: exp.get("n_cal", 3)]:
                for m in methods:
                    if m in exp.get("skip", {}).get(fam_id, ()):
                        continue
                    k = (fam_id, f"{float(p):.4f}", str(seed), m)
                    if k not in done:
                        tasks.append((a.experiment, fam_id, p, seed, m))
    # slow methods first so the pool tail is short
    order = {m: i for i, m in enumerate(exp.get("slow_first", []))}
    tasks.sort(key=lambda t: order.get(t[4], 99))
    print(f"{len(tasks)} tasks to run ({len(done)} already done) -> {path}", flush=True)
    if a.dry:
        return

    new_file = not os.path.exists(path)
    fields = FIELDS
    if not new_file:                      # resume into a file written by an older schema
        with open(path, newline="") as fh:
            fields = next(csv.reader(fh))
    with open(path, "a", newline="") as fh, Pool(a.workers, maxtasksperchild=50) as pool:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new_file:
            w.writeheader()
        t0 = time.time()
        for i, rows in enumerate(pool.imap_unordered(run_task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            if i % 25 == 0 or i == len(tasks):
                el = time.time() - t0
                print(f"  {i}/{len(tasks)}  {el/60:.1f} min  eta {el/i*(len(tasks)-i)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
