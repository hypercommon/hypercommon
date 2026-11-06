"""
Commonality as a family, on the real networks: for every f in the bench family,
how well does f(u,v) on an EDGE separate same-community from different-community
pairs (both endpoints labelled)?

Outputs out/<name>_family.json with, per f:
  auc            pooled AUC, inner vs outer
  auc_strata     AUC within strata of (smaller degree, larger degree) -- controls
                 both endpoint degrees, not just the larger one
  auc_strat_mean edge-count-weighted mean of the stratum AUCs
  box            5/25/50/75/95 percentiles of f for inner and outer edges
and a per-edge sample (out/<name>_family_sample.npz) for figures.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "bench"))

from real_networks import OUT, auc, load_karate, load_snap   # noqa: E402
from hc_family import FUNCTIONS                                # noqa: E402

DEG_EDGES = [1, 3, 5, 8, 12, 20, 35, 60, 100, 10**9]


def edge_table(adj, memb):
    deg = {u: len(v) for u, v in adj.items()}
    rows = []
    for u, nb in adj.items():
        mu = memb.get(u)
        if not mu:
            continue
        for v in nb:
            if u < v and memb.get(v):
                common = adj[u] & adj[v]
                I = len(common) + 2
                su, sv = deg[u] + 1, deg[v] + 1
                ra = sum(1.0 / deg[w] for w in common)
                rows.append((bool(mu & memb[v]), I, su + sv - I, su, sv, ra,
                             min(deg[u], deg[v]), max(deg[u], deg[v])))
    A = np.array(rows, float)
    P = {"I": A[:, 1], "U": A[:, 2], "su": A[:, 3], "sv": A[:, 4], "ra": A[:, 5]}
    return A[:, 0].astype(bool), P, A[:, 6], A[:, 7]


def run(name):
    adj, memb = load_karate() if name == "Karate" else load_snap(name)
    inner, P, dmin, dmax = edge_table(adj, memb)
    bmin = np.digitize(dmin, DEG_EDGES)
    bmax = np.digitize(dmax, DEG_EDGES)
    res = {"name": name, "n_inner": int(inner.sum()), "n_outer": int((~inner).sum()), "f": {}}
    rng = np.random.default_rng(0)
    take = np.arange(len(inner)) if len(inner) <= 200_000 else rng.choice(len(inner), 200_000, replace=False)
    sample = {"inner": inner[take]}
    for fname, f in FUNCTIONS.items():
        x = np.asarray(f(P), float)
        strata, w = [], []
        for a in np.unique(bmin):
            for b in np.unique(bmax[bmin == a]):
                s = (bmin == a) & (bmax == b)
                if (s & inner).sum() >= 20 and (s & ~inner).sum() >= 20:
                    strata.append({"dmin_bin": int(a), "dmax_bin": int(b), "n": int(s.sum()),
                                   "auc": auc(x[s & inner], x[s & ~inner])})
                    w.append(s.sum())
        res["f"][fname] = {
            "auc": auc(x[inner], x[~inner]),
            "auc_strata": strata,
            "auc_strat_mean": float(np.average([st["auc"] for st in strata], weights=w)) if w else None,
            "box": {g: np.percentile(x[m], [5, 25, 50, 75, 95]).tolist()
                    for g, m in (("inner", inner), ("outer", ~inner))},
        }
        sample[fname] = x[take].astype(np.float32)
        print(f"  {name:8s} {fname:15s} AUC {res['f'][fname]['auc']:.3f}  stratified "
              f"{res['f'][fname]['auc_strat_mean'] if w else float('nan'):.3f}", flush=True)
    with open(os.path.join(OUT, f"{name}_family.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    np.savez_compressed(os.path.join(OUT, f"{name}_family_sample.npz"), **sample)


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["Karate", "Dblp", "Amazon", "Youtube"]):
        run(nm)
