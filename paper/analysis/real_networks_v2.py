"""
Second pass on the real networks.

1. Which commonality formula separates same-community edges best, on LABELLED
   pairs only and within degree strata (the thesis compared them on pooled
   classes, where unlabelled nodes reversed the DBLP and YouTube results).

2. Pair rule vs triple rule at matched recall. Wedges are sampled uniformly
   among wedges whose three nodes are all labelled (centre drawn proportional
   to C(labelled degree, 2)), so YouTube -- 4.6% labelled -- is not starved.
   For each wedge we keep the three closed-Jaccard values, so precision/recall
   curves for both rules can be drawn from the same sample.
"""

from __future__ import annotations

import json
import os
import random
import sys

import numpy as np

from real_networks import (OUT, DEG_BINS, auc, load_karate, load_snap,
                           jaccard_closed)

rng = random.Random(11)
WEDGES = 400_000


def edge_measures(adj, memb):
    """Every formula the thesis tried, on edges with both endpoints labelled."""
    deg_mean = np.mean([len(v) for v in adj.values()])
    rows = []
    for u, nb in adj.items():
        mu = memb.get(u)
        if not mu:
            continue
        for v in nb:
            if u < v and memb.get(v):
                inner = bool(mu & memb[v])
                a, b = adj[u], adj[v]
                c = len(a & b) if len(a) < len(b) else len(b & a)
                du, dv = len(a), len(b)
                inter, union = c + 2, du + dv - c          # closed neighbourhoods
                rows.append((inner, inter / union, inter, inter / (min(du, dv) + 1),
                             inter * inter / (deg_mean * union), max(du, dv)))
    A = np.array(rows, float)
    inner = A[:, 0].astype(bool)
    names = ["jaccard", "count", "min_norm", "f3"]
    res = {"n_inner": int(inner.sum()), "n_outer": int((~inner).sum()), "auc": {}, "auc_by_degree": {}}
    for i, nm in enumerate(names, start=1):
        res["auc"][nm] = auc(A[inner, i], A[~inner, i])
        strata = []
        for lo, hi in DEG_BINS:
            s = (A[:, 5] >= lo) & (A[:, 5] < hi)
            if (s & inner).sum() >= 20 and (s & ~inner).sum() >= 20:
                strata.append({"lo": lo, "hi": hi, "n": int(s.sum()),
                               "auc": auc(A[s & inner, i], A[s & ~inner, i])})
        res["auc_by_degree"][nm] = strata
        # degree-weighted average of within-stratum AUC (Mantel-Haenszel style)
        w = np.array([x["n"] for x in strata], float)
        res.setdefault("auc_stratified_mean", {})[nm] = float((w * np.array([x["auc"] for x in strata])).sum() / w.sum()) if len(w) else float("nan")
    return res


def labelled_wedges(adj, memb):
    lab_nb = {u: [v for v in nb if memb.get(v)] for u, nb in adj.items() if memb.get(u)}
    centres = [u for u, l in lab_nb.items() if len(l) >= 2]
    w = np.array([len(lab_nb[u]) * (len(lab_nb[u]) - 1) / 2 for u in centres])
    cum = np.cumsum(w)
    picks = np.searchsorted(cum, np.random.default_rng(5).random(WEDGES) * cum[-1])
    out = np.empty((WEDGES, 5))
    for r, ci in enumerate(picks):
        j = centres[ci]
        i, k = rng.sample(lab_nb[j], 2)
        closed = k in adj[i]
        out[r] = (bool(memb[i] & memb[j] & memb[k]),
                  jaccard_closed(adj, i, j, True),
                  jaccard_closed(adj, j, k, True),
                  jaccard_closed(adj, i, k, closed),
                  closed)
    return out, float(w.sum())


def weighted_auc(pos, neg, wp, wn):
    """AUC with sample weights: P(pos > neg) + 0.5 P(tie)."""
    v = np.concatenate([pos, neg])
    w = np.concatenate([wp, wn])
    lab = np.concatenate([np.ones(len(pos), bool), np.zeros(len(neg), bool)])
    order = np.argsort(v, kind="mergesort")
    v, w, lab = v[order], w[order], lab[order]
    total, cum_neg, i = 0.0, 0.0, 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[j + 1] == v[i]:
            j += 1
        blk_w, blk_l = w[i:j + 1], lab[i:j + 1]
        wpos, wneg = blk_w[blk_l].sum(), blk_w[~blk_l].sum()
        total += wpos * (cum_neg + 0.5 * wneg)
        cum_neg += wneg
        i = j + 1
    return float(total / (wp.sum() * wn.sum()))


def pr_curves(W):
    """
    Precision/recall over UNIQUE triples. Sampling a centre then two of its
    neighbours reaches an open wedge through its one centre but a closed
    triangle through each of its three nodes, so triangles carry weight 1/3.
    """
    pure = W[:, 0].astype(bool)
    closed = W[:, 4].astype(bool)
    wt = np.where(closed, 1.0 / 3.0, 1.0)
    # pair rule = the benchmark's "two-test" ablation: every EDGE of the triple must
    # pass (two for an open wedge, all three for a triangle), so the decision does
    # not depend on which centre the triple was sampled from
    all3 = np.minimum(np.minimum(W[:, 1], W[:, 2]), W[:, 3])
    pair_score = np.where(closed, all3, np.minimum(W[:, 1], W[:, 2]))
    triple_score = all3                                 # all three pairs pass
    grid = np.round(np.arange(0, 1.0001, 0.01), 2)
    curves = {"t": grid.tolist(), "pair": [], "triple": []}
    wpure = (wt * pure).sum()
    for nm, s in (("pair", pair_score), ("triple", triple_score)):
        for t in grid:
            k = s >= t
            wk = wt[k].sum()
            curves[nm].append({"kept": int(k.sum()),
                               "precision": float((wt[k] * pure[k]).sum() / wk) if k.any() else None,
                               "recall": float((wt[k] * pure[k]).sum() / wpure)})
    curves["triangle_only"] = {"precision": float(pure[closed].mean()),
                               "recall": float((wt * (closed & pure)).sum() / wpure)}
    curves["frac_triangles_unique"] = float((wt * closed).sum() / wt.sum())
    curves["base_rate"] = float(wpure / wt.sum())
    curves["auc_pair"] = weighted_auc(pair_score[pure], pair_score[~pure], wt[pure], wt[~pure])
    curves["auc_triple"] = weighted_auc(triple_score[pure], triple_score[~pure], wt[pure], wt[~pure])

    # precision at matched recall R, with a bootstrap over sampled triples
    def prec_at(score, idx, R):
        s, p_, w_ = score[idx], pure[idx], wt[idx]
        order = np.argsort(-s, kind="mergesort")
        cw = np.cumsum((w_ * p_)[order])
        tot = cw[-1]
        k = int(np.searchsorted(cw, R * tot))
        kept = order[:k + 1]
        return float((w_[kept] * p_[kept]).sum() / w_[kept].sum())
    rng_b = np.random.default_rng(3)
    matched = {}
    for R in (0.1, 0.2, 0.3, 0.4):
        full = np.arange(len(W))
        pp, pt = prec_at(pair_score, full, R), prec_at(triple_score, full, R)
        diffs = []
        for _ in range(200):
            b = rng_b.integers(0, len(W), len(W))
            diffs.append(prec_at(triple_score, b, R) - prec_at(pair_score, b, R))
        matched[str(R)] = {"pair": pp, "triple": pt, "diff": pt - pp,
                           "diff_lo": float(np.percentile(diffs, 2.5)), "diff_hi": float(np.percentile(diffs, 97.5))}
    curves["matched_recall"] = matched
    return curves


def run(name):
    print("==", name, flush=True)
    adj, memb = load_karate() if name == "Karate" else load_snap(name)
    res = {"name": name, "edges": edge_measures(adj, memb)}
    print("  measures", {k: round(v, 3) for k, v in res["edges"]["auc"].items()},
          "stratified", {k: round(v, 3) for k, v in res["edges"]["auc_stratified_mean"].items()}, flush=True)
    W, total = labelled_wedges(adj, memb)
    res["labelled_wedges_total"] = total
    res["wedges"] = pr_curves(W)
    print("  wedge auc pair %.3f triple %.3f base %.3f triangle-only P=%.3f R=%.3f" % (
        res["wedges"]["auc_pair"], res["wedges"]["auc_triple"], res["wedges"]["base_rate"],
        res["wedges"]["triangle_only"]["precision"], res["wedges"]["triangle_only"]["recall"]), flush=True)
    np.save(os.path.join(OUT, f"{name}_wedges.npy"), W[:200_000].astype(np.float32))
    with open(os.path.join(OUT, f"{name}_v2.json"), "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["Karate", "Amazon", "Dblp", "Youtube"]):
        run(nm)
