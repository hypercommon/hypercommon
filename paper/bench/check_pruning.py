"""
Degree-ratio pruning is exact: on the full SNAP Amazon and DBLP graphs, the
covers of hc_big.BigScores with pruning at t_min = 0.02 equal those without pruning at
every threshold of the grid, for HyperCommon and its three ablations. Also reports the
wedges enumerated with and without pruning. Writes results/check_pruning.csv.

    python check_pruning.py [amazon dblp]
"""
import csv
import gzip
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hc_big import BigScores  # noqa: E402

SNAP = os.path.join(HERE, "..", "analysis", "snap")
T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]
VARIANTS = ["hc", "hc_tri", "hc_2test", "hc_econn"]


def load(net):
    e = np.loadtxt(gzip.open(os.path.join(SNAP, f"com-{net}.ungraph.txt.gz"), "rt"), comments="#", dtype=np.int64)
    e = e[e[:, 0] != e[:, 1]]
    e = np.unique(np.sort(e, axis=1), axis=0)
    ids = np.unique(e)
    return np.searchsorted(ids, e), len(ids)


def canon(node, comm):
    order = np.lexsort((node, comm))
    node, comm = node[order], comm[order]
    cuts = np.flatnonzero(np.diff(comm)) + 1
    return sorted(tuple(b.tolist()) for b in np.split(node, cuts) if len(b))


def main(nets):
    out = os.path.join(HERE, "results", "check_pruning.csv")
    rows = []
    for net in nets:
        e, n = load(net)
        t0 = time.time()
        full = BigScores(e, n, t_min=0.0)
        pruned = BigScores(e, n, t_min=0.02)
        print(f"{net}: wedges {full.n_wedges} -> {pruned.n_wedges} ({time.time() - t0:.0f} s)", flush=True)
        for v in VARIANTS:
            same = 0
            for t in T_GRID:
                a = canon(*full.cover(t, variant=v))
                b = canon(*pruned.cover(t, variant=v))
                same += a == b
            rows.append(dict(network=net, variant=v, thresholds=len(T_GRID), identical=same,
                             wedges_full=full.n_wedges, wedges_pruned=pruned.n_wedges))
            print(f"  {v}: identical covers at {same}/{len(T_GRID)} thresholds", flush=True)
        del full, pruned
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1:] or ["amazon", "dblp"])
