"""
Real-network evidence for commonality (closed-neighbourhood Jaccard).

Recomputes, from the raw SNAP edge lists and ground-truth community files, the
statistics the thesis plotted -- but with the pair classes defined so that the
comparison is not confounded by unlabelled nodes:

  inner      both endpoints labelled, share >= 1 ground-truth community
  outer      both endpoints labelled, share none
  unlabelled at least one endpoint belongs to no ground-truth community

SNAP ground truth covers only part of each network, and an unlabelled node is
"no information", not "different community". Pooling it into outer is what made
DBLP look inverted in the thesis.

Outputs JSON summaries (and per-pair samples for plotting) under
paper/analysis/out/.
"""

from __future__ import annotations

import json
import math
import os
import random
import sys
from collections import defaultdict

import networkx as nx
import numpy as np
from scipy.stats import rankdata

# Public SNAP files (https://snap.stanford.edu/data/#communities), downloaded to
# paper/analysis/snap/. Verified identical to the copies used in the thesis.
SNAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "snap")
SNAP_FILES = {  # name -> (graph, all ground-truth communities)
    "Dblp": ("com-dblp.ungraph.txt.gz", "com-dblp.all.cmty.txt.gz"),
    "Amazon": ("com-amazon.ungraph.txt.gz", "com-amazon.all.dedup.cmty.txt.gz"),
    "Youtube": ("com-youtube.ungraph.txt.gz", "com-youtube.all.cmty.txt.gz"),
}
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)

T_GRID = [round(0.05 * i, 2) for i in range(0, 21)]
DEG_BINS = [(1, 5), (5, 10), (10, 20), (20, 50), (50, 100), (100, 10**9)]
WEDGE_SAMPLES = 2_000_000
rng = random.Random(20261022)


# ---------------------------------------------------------------- loading

def load_snap(name):
    import gzip
    graph_file, cmty_file = SNAP_FILES[name]
    adj = defaultdict(set)
    with gzip.open(os.path.join(SNAP, graph_file), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, b = line.split()[:2]
            a, b = int(a), int(b)
            if a != b:
                adj[a].add(b)
                adj[b].add(a)
    memb = defaultdict(set)
    with gzip.open(os.path.join(SNAP, cmty_file), "rt") as fh:
        for cid, line in enumerate(fh):
            for x in line.split():
                memb[int(x)].add(cid)
    return dict(adj), dict(memb)


def load_karate():
    G = nx.karate_club_graph()
    adj = {u: set(G.neighbors(u)) for u in G}
    memb = {u: {0 if G.nodes[u]["club"] == "Mr. Hi" else 1} for u in G}
    return adj, memb


# ---------------------------------------------------------------- measures

def jaccard_closed(adj, u, v, adjacent):
    """|N[u] & N[v]| / |N[u] | N[v]| with closed neighbourhoods."""
    a, b = adj[u], adj[v]
    c = len(a & b) if len(a) < len(b) else len(b & a)
    du, dv = len(a), len(b)
    if adjacent:
        # u in N[v] and v in N[u]: both endpoints are in the intersection
        return (c + 2) / (du + dv - c)
    return c / (du + dv + 2 - c)


def pair_class(memb, u, v):
    mu, mv = memb.get(u), memb.get(v)
    if not mu or not mv:
        return "unlabelled"
    return "inner" if mu & mv else "outer"


def auc(pos, neg):
    """P(score_pos > score_neg) + 0.5 P(tie), via the rank-sum statistic."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = rankdata(np.concatenate([pos, neg]))  # average ranks for ties
    r_pos = ranks[: len(pos)].sum()
    return float((r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def quartiles(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return None
    q = np.percentile(x, [5, 25, 50, 75, 95])
    return {"n": int(len(x)), "mean": float(x.mean()),
            "p05": q[0], "q1": q[1], "median": q[2], "q3": q[3], "p95": q[4]}


# ---------------------------------------------------------------- analyses

def edge_analysis(adj, memb):
    cls, js, dmax, dmin, tri = [], [], [], [], []
    for u, nb in adj.items():
        for v in nb:
            if u < v:
                c = pair_class(memb, u, v)
                cls.append(c)
                js.append(jaccard_closed(adj, u, v, True))
                du, dv = len(adj[u]), len(adj[v])
                dmax.append(max(du, dv))
                dmin.append(min(du, dv))
                tri.append(len(adj[u] & adj[v]))
    cls = np.array(cls)
    js = np.array(js)
    dmax = np.array(dmax)
    tri = np.array(tri)

    res = {"edges": int(len(cls)),
           "class_counts": {k: int((cls == k).sum()) for k in ("inner", "outer", "unlabelled")}}
    res["jaccard"] = {k: quartiles(js[cls == k]) for k in ("inner", "outer", "unlabelled")}
    res["jaccard"]["outer_pooled"] = quartiles(js[cls != "inner"])
    res["auc_inner_vs_outer"] = auc(js[cls == "inner"], js[cls == "outer"])
    res["auc_inner_vs_pooled"] = auc(js[cls == "inner"], js[cls != "inner"])

    # Degree-stratified: does inner > outer hold at fixed degree?
    strat = []
    for lo, hi in DEG_BINS:
        s = (dmax >= lo) & (dmax < hi)
        si, so = s & (cls == "inner"), s & (cls == "outer")
        if si.sum() < 20 or so.sum() < 20:
            continue
        strat.append({"dmax_lo": lo, "dmax_hi": hi,
                      "n_inner": int(si.sum()), "n_outer": int(so.sum()),
                      "median_inner": float(np.median(js[si])),
                      "median_outer": float(np.median(js[so])),
                      "auc": auc(js[si], js[so])})
    res["degree_stratified"] = strat

    # Triangle support: edges with no common neighbour can never be in a triple
    lab = cls != "unlabelled"
    res["frac_edges_no_triangle"] = {k: float((tri[cls == k] == 0).mean()) for k in ("inner", "outer")}

    # Pair-rule precision / recall over the threshold grid (labelled edges only)
    pr = []
    n_inner = (cls == "inner").sum()
    for t in T_GRID:
        keep = lab & (js >= t)
        tp = (keep & (cls == "inner")).sum()
        fp = (keep & (cls == "outer")).sum()
        pr.append({"t": t, "kept": int(keep.sum()),
                   "precision": float(tp / (tp + fp)) if tp + fp else float("nan"),
                   "recall": float(tp / n_inner) if n_inner else float("nan")})
    res["pair_rule"] = pr
    res["base_rate_inner"] = float(n_inner / lab.sum())

    # small sample for plotting
    idx = np.arange(len(cls))
    if len(idx) > 300_000:
        idx = np.array(rng.sample(range(len(cls)), 300_000))
    sample = {"j": js[idx].round(5).tolist(), "cls": cls[idx].tolist(), "dmax": dmax[idx].tolist()}
    return res, sample


def wedge_analysis(adj, memb, exact=False):
    """
    Triples with a centre (a wedge i-j-k, i and k both neighbours of j) --
    the objects HyperCommon admits. A triple is `pure` when one ground-truth
    community contains all three nodes, `mixed` when all three are labelled
    but no community holds them all.

    Compares the triple rule (all three pairs >= t) against the pair rule on
    the two edges of the wedge, to measure what requiring a triangle of
    agreement buys in precision.
    """
    nodes = [u for u in adj if len(adj[u]) >= 2]
    weights = np.array([len(adj[u]) * (len(adj[u]) - 1) / 2 for u in nodes], float)
    total_wedges = float(weights.sum())

    def gen():
        if exact:
            for j in nodes:
                nb = sorted(adj[j])
                for a in range(len(nb)):
                    for b in range(a + 1, len(nb)):
                        yield nb[a], j, nb[b]
        else:
            cum = np.cumsum(weights)
            centres = np.searchsorted(cum, np.random.default_rng(7).random(WEDGE_SAMPLES) * cum[-1])
            for ci in centres:
                j = nodes[ci]
                i, k = rng.sample(tuple(adj[j]), 2)
                yield i, j, k

    rows = []
    for i, j, k in gen():
        mi, mj, mk = memb.get(i), memb.get(j), memb.get(k)
        if not mi or not mj or not mk:
            continue
        pure = bool(mi & mj & mk)
        ik_adj = k in adj[i]
        jij = jaccard_closed(adj, i, j, True)
        jjk = jaccard_closed(adj, j, k, True)
        jik = jaccard_closed(adj, i, k, ik_adj)
        rows.append((pure, jij, jjk, jik, ik_adj))

    arr = np.array([(r[0], r[1], r[2], r[3], r[4]) for r in rows], float)
    pure = arr[:, 0].astype(bool)
    jmin_edges = np.minimum(arr[:, 1], arr[:, 2])
    jmin_all = np.minimum(jmin_edges, arr[:, 3])
    closed = arr[:, 4].astype(bool)

    out = {"total_wedges": total_wedges, "evaluated_labelled": int(len(arr)),
           "base_rate_pure": float(pure.mean()), "sampled": not exact,
           "frac_closed_pure": float(closed[pure].mean()) if pure.any() else None,
           "frac_closed_mixed": float(closed[~pure].mean()) if (~pure).any() else None,
           "auc_min_all": auc(jmin_all[pure], jmin_all[~pure]),
           "auc_min_edges": auc(jmin_edges[pure], jmin_edges[~pure])}
    curve = []
    for t in T_GRID:
        a = jmin_edges >= t          # both wedge edges pass: what a pair rule allows
        b = jmin_all >= t            # all three pairs pass: the triple rule
        curve.append({"t": t,
                      "pair_rule_kept": int(a.sum()),
                      "pair_rule_precision": float(pure[a].mean()) if a.any() else float("nan"),
                      "triple_rule_kept": int(b.sum()),
                      "triple_rule_precision": float(pure[b].mean()) if b.any() else float("nan"),
                      "triple_rule_recall": float((b & pure).sum() / pure.sum()) if pure.any() else float("nan")})
    out["curve"] = curve
    return out


def basic_stats(adj, memb):
    deg = np.array([len(v) for v in adj.values()])
    labelled = sum(1 for u in adj if memb.get(u))
    multi = sum(1 for u in adj if len(memb.get(u, ())) > 1)
    return {"nodes": int(len(adj)), "edges": int(deg.sum() // 2),
            "mean_degree": float(deg.mean()), "median_degree": float(np.median(deg)),
            "max_degree": int(deg.max()),
            "frac_labelled": labelled / len(adj),
            "frac_multi_membership_of_labelled": multi / labelled if labelled else 0.0,
            "communities": len({c for s in memb.values() for c in s})}


def run(name):
    print(f"== {name}", flush=True)
    adj, memb = load_karate() if name == "Karate" else load_snap(name)
    res = {"name": name, "basic": basic_stats(adj, memb)}
    print("  basic", res["basic"], flush=True)
    res["edges"], sample = edge_analysis(adj, memb)
    print("  edges done", flush=True)
    res["wedges"] = wedge_analysis(adj, memb, exact=(name == "Karate"))
    print("  wedges done", flush=True)
    with open(os.path.join(OUT, f"{name}.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    with open(os.path.join(OUT, f"{name}_edge_sample.json"), "w") as fh:
        json.dump(sample, fh)


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["Karate", "Amazon", "Dblp", "Youtube"]):
        run(nm)
