"""
The nested degree-corrected stochastic block model (Peixoto 2014; graph-tool
minimize_nested_blockmodel_dl) on the benchmark graphs, a baseline that selects its own
number of groups by minimum description length.

Runs under Linux (WSL) with graph-tool, on the same graphs as run2.py (bench/graphs.py,
string-seeded; LFR graphs from bench/lfr_cache). Every level of the fitted hierarchy is
a partition and is scored like any other cover, with the level index as the parameter
(0 = finest), so that the protocols read it as for the other methods: oracle = best
level, label-free = level of largest overlapping modularity, native = level 0 (the
model's own finest partition).

    python run_sbm.py <experiment> [--workers N]
Rows go to results2/<experiment>_sbm.csv (method "sbm"), appended and resumable.
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
from run2 import CAL_SEEDS, FIELDS, TEST_SEEDS, _safe  # noqa: E402

RESULTS = os.path.join(HERE, "results2")


def run_task(args):
    exp_name, fam_id, p, seed = args
    import numpy as np
    import graph_tool.all as gt
    import cover_metrics as M
    from self_tune import _adj, overlapping_modularity
    exp = EXPERIMENTS[exp_name]
    fam = exp["families"][fam_id]
    split = "cal" if seed >= 1000 else "test"
    rows, t0 = [], time.time()
    try:
        G, truth = exp["build"](fam, p, seed)
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        g = gt.Graph(directed=False)
        g.add_vertex(len(nodes))
        g.add_edge_list([(idx[u], idx[v]) for u, v in G.edges()])
        gt.seed_rng(seed % (2**31) + 1)
        np.random.seed(seed % (2**31) + 1)
        state = gt.minimize_nested_blockmodel_dl(g)
        A, aidx = _adj(G)
        sizes = np.array([len(c) for c in truth])
        med = np.median(sizes)
        small, large = sizes <= med, sizes > med
        fit_s = time.time() - t0
        for level in range(len(state.levels)):
            b = state.project_level(level).get_blocks().a
            groups = {}
            for i, lab in enumerate(b):
                groups.setdefault(int(lab), set()).add(nodes[i])
            cover = list(groups.values())
            if len(cover) <= 1 and level > 0:
                break
            rec = M.truth_recovery(cover, truth, G.nodes())
            rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split, method="sbm",
                             param=level, omega=M.omega(cover, truth, G.nodes()),
                             onmi=_safe(M.onmi, cover, truth, G), f1=M.f1_best_match(cover, truth, G.nodes()),
                             coverage=M.coverage(cover, truth), n_comms=len(cover),
                             eq=overlapping_modularity(A, aidx, cover), rec_mean=float(rec.mean()),
                             rec_small=float(rec[small].mean()) if small.any() else "",
                             rec_large=float(rec[large].mean()) if large.any() else "",
                             n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(), n_wedges="",
                             seconds=round(fit_s, 3), error=""))
    except Exception as e:
        rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split, method="sbm",
                         seconds=round(time.time() - t0, 3),
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    exp = EXPERIMENTS[a.experiment]
    path = os.path.join(RESULTS, f"{a.experiment}_sbm.csv")
    done = set()
    if os.path.exists(path):
        with open(path, newline="") as fh:
            done = {(r["family"], f"{float(r['p']):.4f}", str(r["seed"])) for r in csv.DictReader(fh)}
    tasks = [(a.experiment, f, p, s) for f in exp["families"]
             for p in exp.get("p_grid_by_family", {}).get(f, exp.get("p_grid_v2", exp["p_grid"]))
             for s in TEST_SEEDS + CAL_SEEDS if (f, f"{float(p):.4f}", str(s)) not in done]
    print(f"{len(tasks)} tasks ({len(done)} done) -> {path}", flush=True)
    if a.dry:
        return
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh, Pool(a.workers, maxtasksperchild=20) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        t0 = time.time()
        for i, rows in enumerate(pool.imap_unordered(run_task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            if i % 50 == 0 or i == len(tasks):
                el = time.time() - t0
                print(f"  {i}/{len(tasks)}  {el / 60:.1f} min  eta {el / i * (len(tasks) - i) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
