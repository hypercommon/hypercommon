"""
HyperCommon with a FAMILY of commonality functions.

The method only needs a symmetric pair score f(u,v); Jaccard is one choice. For
every wedge we keep the three pairs it involves, and for every pair the raw
quantities any neighbourhood-overlap score is built from:

  I   = |N[u] & N[v]|       closed neighbourhoods (u, v count when adjacent)
  U   = |N[u] | N[v]|
  su  = |N[u]|, sv = |N[v]|
  ra  = sum over common open neighbours w of 1/deg(w)   (resource allocation)

A triple's score under f is the minimum of f over its three pairs (all three must
pass), and communities at a threshold are the components of the shared-pair
graph on admitted triples, exactly as for Jaccard.

Different f live on different scales, so thresholds are taken as QUANTILES q of
the f-scores of all candidate triples of the graph: the same q grid gives every
f the same number of operating points, spread over its own range.
"""

from __future__ import annotations

import numpy as np

from hc_fast import TripleScores

FUNCTIONS = {
    "jaccard": lambda P: P["I"] / P["U"],
    "count": lambda P: P["I"].astype(float),
    "i2_over_u": lambda P: P["I"] ** 2 / P["U"],          # thesis f3 up to a constant
    "i_over_sqrtu": lambda P: P["I"] / np.sqrt(P["U"]),
    "cosine": lambda P: P["I"] / np.sqrt(P["su"] * P["sv"]),
    "min_norm": lambda P: P["I"] / np.minimum(P["su"], P["sv"]),
    "max_norm": lambda P: P["I"] / np.maximum(P["su"], P["sv"]),
    "resource_alloc": lambda P: P["ra"],
}
Q_GRID = [round(0.02 * i, 2) for i in range(0, 50)]        # quantile levels 0 .. 0.98


class FamilyScores(TripleScores):
    def __init__(self, G):
        super().__init__(G)
        nodes = self.nodes
        idx = self.index
        nbr = [set(idx[v] for v in G.neighbors(u)) for u in nodes]
        deg = np.array([len(s) for s in nbr], float)
        t = self.triples
        pairs = np.concatenate([t[:, [0, 1]], t[:, [0, 2]], t[:, [1, 2]]])
        pkey = pairs[:, 0] * self.n + pairs[:, 1]
        uniq = np.unique(pkey)
        a_ids, b_ids = uniq // self.n, uniq % self.n
        I = np.empty(len(uniq), np.int64)
        ra = np.empty(len(uniq))
        adj = np.empty(len(uniq), bool)
        for r, (a, b) in enumerate(zip(a_ids, b_ids)):
            common = nbr[a] & nbr[b]
            adj[r] = b in nbr[a]
            I[r] = len(common) + (2 if adj[r] else 0)
            ra[r] = sum(1.0 / deg[w] for w in common)
        su, sv = deg[a_ids] + 1, deg[b_ids] + 1
        self.P = {"I": I, "U": su + sv - I, "su": su, "sv": sv, "ra": ra}
        self._fcache = {}

    def triple_score(self, fname):
        s = self._fcache.get(fname)
        if s is None:
            pv = FUNCTIONS[fname](self.P)
            s = pv[self.pair_of].min(axis=1)        # all three pairs must pass
            self._fcache[fname] = s
        return s

    def communities_f(self, fname, q):
        s = self.triple_score(fname)
        thr = np.quantile(s, q) if q > 0 else -np.inf
        # reuse the component machinery through the s3 slot
        saved = self.s3
        self.s3 = s
        try:
            return self.communities("hc", thr)
        finally:
            self.s3 = saved
