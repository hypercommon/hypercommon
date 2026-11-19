"""
How much planted structure is left at each rewiring level p, on the ring families of the
robustness sweep (no detection, graphs only). Per family, p in 0, 0.01, ..., 1 and test
seeds 0..9:

  inside       mean over planted rings of the fraction of their members' edge ends that
               stay inside the ring, in/(in+out)  (about 1 - p under uniform rewiring)
  weak         fraction of planted rings that are communities in the weak sense of
               Radicchi et al. (2004): more edge ends inside than outside, in > out
  strong       fraction of nodes with more neighbours inside one of their rings than
               outside all of them
  clustering   average clustering coefficient, and its ratio to the same graph at p = 0
               (the lattice signal HyperCommon relies on; about (1-p)^3 for one ring,
               Barrat and Weigt 2000)

    python robust_structure.py [--workers N]
Writes results2/robust_structure.csv (one row per family, p, seed).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from experiments import EXPERIMENTS  # noqa: E402
from run2 import P_DENSE, TEST_SEEDS  # noqa: E402

EXPS = ["rings_uniform", "rings_hetero", "rings_hetero_overlap", "rings_pref"]
OUT = os.path.join(HERE, "results2", "robust_structure.csv")
FIELDS = ["experiment", "family", "p", "seed", "inside", "weak", "strong", "clustering"]


def task(args):
    import networkx as nx
    e, fam_id, p, seed = args
    G, truth = EXPERIMENTS[e]["build"](EXPERIMENTS[e]["families"][fam_id], p, seed)
    member = {}
    for c, com in enumerate(truth):
        for u in com:
            member.setdefault(u, set()).add(c)
    ins, wk = [], 0
    for c, com in enumerate(truth):
        i = o = 0
        for u in com:
            for v in G.adj[u]:
                if v in com:
                    i += 1
                else:
                    o += 1
        ins.append(i / (i + o) if i + o else 0.0)
        wk += i > o
    strong = 0
    for u, cs in member.items():
        nb = G.adj[u]
        best = max(sum(1 for v in nb if c in member.get(v, ())) for c in cs)
        strong += best > len(nb) - best
    return dict(experiment=e, family=fam_id, p=p, seed=seed, inside=sum(ins) / len(ins),
                weak=wk / len(truth), strong=strong / len(member),
                clustering=nx.average_clustering(G))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    tasks = [(e, f, p, s) for e in EXPS for f in EXPERIMENTS[e]["families"] for p in P_DENSE for s in TEST_SEEDS]
    if a.limit:
        tasks = tasks[:a.limit]
    out = OUT if not a.limit else OUT.replace(".csv", "_test.csv")
    print(f"{len(tasks)} graphs -> {out}", flush=True)
    with open(out, "w", newline="") as fh, Pool(a.workers) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for i, row in enumerate(pool.imap_unordered(task, tasks, chunksize=4), 1):
            w.writerow(row)
            if i % 1000 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)


if __name__ == "__main__":
    main()
