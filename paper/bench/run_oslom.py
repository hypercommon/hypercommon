"""
OSLOM (Lancichinetti, Radicchi, Ramasco & Fortunato 2011; OSLOM2, oslom_undir) on the
benchmark graphs: an overlapping method that keeps a group only if it is statistically
significant against a configuration-model null.

Runs with the compiled OSLOM2 binary (bench/build_oslom.sh; Linux or macOS), on the
same graphs as run2.py (bench/graphs.py, string-seeded; LFR graphs from bench/lfr_cache).
The parameter is OSLOM's significance threshold -t (P-value; default 0.1), swept over
T_OSLOM; every other option is OSLOM's default (-r 10 runs, -cp 0.5) except -hr 0: only
the first hierarchical level is kept, the one that OSLOM reports as its cover. The
protocols read it as for the other methods: oracle = best -t, label-free = -t of largest
overlapping modularity, native = OSLOM's default -t 0.1.

    python run_oslom.py <experiment> [--workers N]
Rows go to results2/<experiment>_oslom.csv (method "oslom"), appended and resumable.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from experiments import EXPERIMENTS  # noqa: E402
from run2 import CAL_SEEDS, FIELDS, TEST_SEEDS, _safe  # noqa: E402

RESULTS = os.path.join(HERE, "results2")
OSLOM_BIN = os.environ.get("OSLOM_BIN", os.path.expanduser("~/tools/OSLOM2/oslom_undir"))
T_OSLOM = [1e-14, 1e-11, 1e-8, 1e-5, 1e-2, 0.1]


def oslom_cover_edges(edges, t, seed, workdir, timeout=None):
    """One OSLOM run on an edge list of integer ids at threshold t; returns (cover as list
    of sets of those ids, seconds). Past `timeout` seconds subprocess.TimeoutExpired."""
    net = os.path.join(workdir, f"net_{t}.dat")
    with open(net, "w") as fh:
        for u, v in edges:
            fh.write(f"{u}\t{v}\n")
    t0 = time.time()
    # relative name (cwd=workdir): OSLOM keeps the file name in a fixed-size buffer, and
    # the long macOS temporary path overflows it (SIGABRT)
    subprocess.run([OSLOM_BIN, "-f", os.path.basename(net), "-uw", "-r", "10", "-hr", "0", "-t", str(t),
                    "-seed", str(seed % (2**31) + 1)],
                   cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
                   timeout=timeout)
    secs = time.time() - t0
    cover = []
    with open(os.path.join(net + "_oslo_files", "tp")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cover.append({int(x) for x in line.split()})
    return cover, secs


def oslom_cover(G, t, seed, workdir):
    """One OSLOM run on G at threshold t; returns (cover as list of node sets, seconds)."""
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    cover, secs = oslom_cover_edges(((idx[u], idx[v]) for u, v in G.edges()), t, seed, workdir)
    return [{nodes[i] for i in c} for c in cover], secs


def run_task(args):
    exp_name, fam_id, p, seed = args
    import numpy as np
    import cover_metrics as M
    from self_tune import _adj, overlapping_modularity
    exp = EXPERIMENTS[exp_name]
    fam = exp["families"][fam_id]
    split = "cal" if seed >= 1000 else "test"
    rows, t0 = [], time.time()
    workdir = tempfile.mkdtemp(prefix="oslom_")
    try:
        G, truth = exp["build"](fam, p, seed)
        A, aidx = _adj(G)
        sizes = np.array([len(c) for c in truth])
        med = np.median(sizes)
        small, large = sizes <= med, sizes > med
        for t in T_OSLOM:
            cover, secs = oslom_cover(G, t, seed, workdir)
            rec = M.truth_recovery(cover, truth, G.nodes())
            rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split, method="oslom",
                             param=t, omega=M.omega(cover, truth, G.nodes()),
                             onmi=_safe(M.onmi, cover, truth, G), f1=M.f1_best_match(cover, truth, G.nodes()),
                             coverage=M.coverage(cover, truth), n_comms=len(cover),
                             eq=overlapping_modularity(A, aidx, cover), rec_mean=float(rec.mean()),
                             rec_small=float(rec[small].mean()) if small.any() else "",
                             rec_large=float(rec[large].mean()) if large.any() else "",
                             n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(), n_wedges="",
                             seconds=round(secs, 3), error=""))
    except Exception as e:
        rows.append(dict(experiment=exp_name, family=fam_id, p=p, seed=seed, split=split, method="oslom",
                         seconds=round(time.time() - t0, 3),
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=2)!r}"[:500]))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="run only the first N tasks (timing)")
    a = ap.parse_args()
    exp = EXPERIMENTS[a.experiment]
    path = os.path.join(RESULTS, f"{a.experiment}_oslom.csv")
    done = set()
    if os.path.exists(path):
        with open(path, newline="") as fh:
            done = {(r["family"], f"{float(r['p']):.4f}", str(r["seed"])) for r in csv.DictReader(fh)}
    tasks = [(a.experiment, f, p, s) for f in exp["families"]
             for p in exp.get("p_grid_by_family", {}).get(f, exp.get("p_grid_v2", exp["p_grid"]))
             for s in TEST_SEEDS + CAL_SEEDS if (f, f"{float(p):.4f}", str(s)) not in done]
    if a.limit:
        tasks = tasks[:a.limit]
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
