"""
Second benchmark round for the commonality family: HyperCommon with each function of
hc_family.FUNCTIONS used in the paper, on the ring and LFR families of the main
tables, every cover scored as in run2.py (ONMI omitted) so that analysis2.py gives the
oracle, frozen and label-free choice of the quantile level for every function.

Rings are sampled every 0.1 in p (the main tables use 0.05), the quantile grid every
0.04 (hc_family.Q_GRID uses 0.02); the Jaccard rows give the same covers as `hc` at
the corresponding thresholds and serve as a check.

    python run2_family.py [--workers N] [--dry] [--top]
Rows go to results2/family2.csv (appended; a rerun skips tasks already present).
With --top only the levels Q_TOP are run, into results2/family2_top.csv: at strong
noise the best level of the main grid is its last one, 0.96.
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

from experiments import EXPERIMENTS, MU_GRID, _lfr, _rings  # noqa: E402
from run2 import CAL_SEEDS, FIELDS, RESULTS, TEST_SEEDS      # noqa: E402

EXP = "family2"
FUNCS = ("jaccard", "cosine", "max_norm", "min_norm", "i2_over_u", "count", "resource_alloc")
Q_GRID = [round(0.04 * i, 2) for i in range(25)]             # 0 .. 0.96
Q_TOP = [0.97, 0.98, 0.99, 0.995]                              # --top: the high end, results2/family2_top.csv
RING_P = [round(0.1 * i, 1) for i in range(9)]               # 0 .. 0.8
FAMILIES = {}
for src, kind in (("rings_uniform", "rings"), ("rings_hetero", "rings"), ("lfr", "lfr"), ("lfr_slices", "lfr")):
    for fam_id, fam in EXPERIMENTS[src]["families"].items():
        FAMILIES[fam_id] = dict(fam, kind=kind, source=src)


def build(fam, p, seed):
    return (_rings if fam["kind"] == "rings" else _lfr)(fam, p, seed)


def _omega_fn(G, truth):
    """Omega against the truth: score_big's grouped computation (equal to
    cover_metrics.omega, test_score_big.py) when the truth covers every node."""
    import numpy as np
    import scipy.sparse as sp
    import cover_metrics as M
    import score_big as SB
    nodes = list(G.nodes())
    index = {u: i for i, u in enumerate(nodes)}
    T = SB.Truth(None, comms=[np.array(sorted(index[u] for u in c)) for c in truth])
    if T.U != len(nodes):
        return lambda cover: M.omega(cover, truth, nodes)
    prep = T.omega_prep["raw"]

    def om(cover):
        cover = [c for c in cover if c]
        if not cover:
            return SB.omega(None, prep)
        r = np.fromiter((index[u] for c in cover for u in c), np.int64)
        c = np.repeat(np.arange(len(cover)), [len(x) for x in cover])
        return SB.omega(sp.csr_matrix((np.ones(len(r), np.int64), (r, c)), shape=(T.U, len(cover))), prep)
    return om


def run_task(args):
    """All functions on one graph (the pair quantities are computed once)."""
    fam_id, p, seed, funcs, qs = args
    import numpy as np
    import cover_metrics as M
    from hc_family import FamilyScores
    from self_tune import _adj, overlapping_modularity

    fam = FAMILIES[fam_id]
    split = "cal" if seed >= 1000 else "test"
    rows = []
    t0 = time.time()
    try:
        G, truth = build(fam, p, seed)
        fs = FamilyScores(G)
        A, idx = _adj(G)
        omega = _omega_fn(G, truth)
        sizes = np.array([len(c) for c in truth])
        med = np.median(sizes)
        small, large = sizes <= med, sizes > med
    except Exception as e:  # recorded, never swallowed silently
        return [dict(experiment=EXP, family=fam_id, p=p, seed=seed, split=split, method="hcf_" + f,
                     seconds=round(time.time() - t0, 3),
                     error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]) for f in funcs]
    t_build = time.time() - t0
    for f in funcs:
        rows += _run_f(fam_id, p, seed, split, f, G, truth, fs, A, idx, omega, small, large, t_build, qs)
        t_build = 0.0
    return rows


def _run_f(fam_id, p, seed, split, f, G, truth, fs, A, idx, omega, small, large, t_build, qs):
    import cover_metrics as M
    from self_tune import overlapping_modularity
    method = "hcf_" + f
    rows = []
    t0 = time.time()
    try:
        t1 = time.time()
        for q in qs:
            cover = fs.communities_f(f, q)
            rec = M.truth_recovery(cover, truth, G.nodes())
            rows.append(dict(experiment=EXP, family=fam_id, p=p, seed=seed, split=split, method=method, param=q,
                             omega=omega(cover), onmi="",
                             f1=M.f1_best_match(cover, truth, G.nodes()),
                             coverage=M.coverage(cover, truth), n_comms=len(cover),
                             eq=overlapping_modularity(A, idx, cover),
                             rec_mean=float(rec.mean()),
                             rec_small=float(rec[small].mean()) if small.any() else "",
                             rec_large=float(rec[large].mean()) if large.any() else "",
                             n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(), n_wedges=fs.n_wedges,
                             seconds=round(time.time() - t1 + t_build, 3), error=""))
            t1, t_build = time.time(), 0.0
    except Exception as e:  # recorded, never swallowed silently
        rows.append(dict(experiment=EXP, family=fam_id, p=p, seed=seed, split=split, method=method,
                         seconds=round(time.time() - t0, 3),
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--families", default="")
    ap.add_argument("--top", action="store_true")
    a = ap.parse_args()
    qs = Q_TOP if a.top else Q_GRID
    path = os.path.join(RESULTS, f"{EXP}_top.csv" if a.top else f"{EXP}.csv")
    done = set()
    if os.path.exists(path):
        with open(path, newline="") as fh:
            done = {(r["family"], f"{float(r['p']):.4f}", str(r["seed"]), r["method"]) for r in csv.DictReader(fh)}
    fams = a.families.split(",") if a.families else list(FAMILIES)
    tasks = []
    for fam_id in fams:
        for p in (RING_P if FAMILIES[fam_id]["kind"] == "rings" else MU_GRID):
            for seed in TEST_SEEDS + CAL_SEEDS:
                todo = tuple(f for f in FUNCS if (fam_id, f"{float(p):.4f}", str(seed), "hcf_" + f) not in done)
                if todo:
                    tasks.append((fam_id, p, seed, todo, qs))
    print(f"{len(tasks)} tasks to run ({len(done)} already done) -> {path}", flush=True)
    if a.dry:
        return
    new_file = not os.path.exists(path)
    with open(path, "a", newline="") as fh, Pool(a.workers, maxtasksperchild=20) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new_file:
            w.writeheader()
        t0 = time.time()
        for i, rows in enumerate(pool.imap_unordered(run_task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            if i % 10 == 0 or i == len(tasks):
                el = time.time() - t0
                print(f"  {i}/{len(tasks)}  {el/60:.1f} min  eta {el/i*(len(tasks)-i)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
