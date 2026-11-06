"""
Threshold sweeps of HyperCommon and its two structural ablations from one pass.

Every wedge (three nodes, one adjacent to the other two) is enumerated once and
given three scores computed from closed-neighbourhood Jaccard J:

  s3      min J over all three pairs               -> HyperCommon (the method)
  s2      min J over the wedge's edges only        -> "two-test" ablation: the
          (for a triangle, all three are edges)       third, non-adjacent pair
                                                      of an open wedge is ignored
  closed  whether the triple is a triangle         -> "triangles only" ablation
                                                      (= 3-clique percolation on
                                                      the graph with edges J < t
                                                      deleted, J from the
                                                      original graph)

A variant admits the triples whose score is >= t (and, for triangles-only, that
are closed). Communities are the node unions of the connected components of the
graph on admitted triples where two triples are adjacent iff they share a pair.
That is computed as components of the bipartite triple--pair incidence graph,
which has the same triple components (two triples sharing a pair are joined
through that pair's node).

`communities(t)` for variant "hc" equals hypercommon.algorithm.get_communities
with closed_neighborhood_jaccard_predicate(t); test_hc_fast.py checks this.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components


class TripleScores:
    def __init__(self, G):
        nodes = list(G.nodes())
        self.index = {u: i for i, u in enumerate(nodes)}
        self.nodes = nodes
        nbr = [set(self.index[v] for v in G.neighbors(u)) for u in nodes]
        deg = np.array([len(s) for s in nbr])
        jcache: dict[tuple[int, int], float] = {}

        def J(a, b):
            if a > b:
                a, b = b, a
            key = (a, b)
            v = jcache.get(key)
            if v is None:
                na, nb = nbr[a], nbr[b]
                c = len(na & nb) if len(na) < len(nb) else len(nb & na)
                if b in na:          # adjacent: both endpoints in the intersection
                    v = (c + 2) / (deg[a] + deg[b] - c)
                else:
                    v = c / (deg[a] + deg[b] + 2 - c)
                jcache[key] = v
            return v

        tri, s3, s2, closed, jp, centre = [], [], [], [], [], []
        seen = set()
        for j in range(len(nodes)):
            nb = sorted(nbr[j])
            for x in range(len(nb)):
                i = nb[x]
                for y in range(x + 1, len(nb)):
                    k = nb[y]
                    a, b, c = sorted((i, j, k))
                    key = (a, b, c)
                    if key in seen:
                        continue
                    is_closed = k in nbr[i]
                    if is_closed:
                        seen.add(key)       # reachable from all three centres
                    jij, jjk, jik = J(i, j), J(j, k), J(i, k)
                    tri.append(key)
                    s3.append(min(jij, jjk, jik))
                    s2.append(min(jij, jjk, jik) if is_closed else min(jij, jjk))
                    closed.append(is_closed)
                    jp.append((J(a, b), J(a, c), J(b, c)))
                    centre.append(-1 if is_closed else j)
        self.triples = np.array(tri, dtype=np.int64).reshape(-1, 3)
        self.s3 = np.array(s3)
        self.s2 = np.array(s2)
        self.closed = np.array(closed, dtype=bool)
        self.jpair = np.array(jp, dtype=float).reshape(-1, 3)     # J of pairs ab, ac, bc
        self.centre = np.array(centre, dtype=np.int64)            # open wedge: its centre
        self.n = len(nodes)
        self.n_wedges = int(sum(d * (d - 1) // 2 for d in deg))
        self.m_edges = int(deg.sum() // 2)

        # pair ids for the incidence graph, built once
        t = self.triples
        pairs = np.concatenate([t[:, [0, 1]], t[:, [0, 2]], t[:, [1, 2]]])
        pkey = pairs[:, 0] * self.n + pairs[:, 1]
        uniq, inv = np.unique(pkey, return_inverse=True)
        self.n_pairs = len(uniq)
        m = len(t)
        self.pair_of = inv.reshape(3, m).T          # (m, 3) pair ids per triple (ab, ac, bc)
        self.pair_nodes = np.stack([uniq // self.n, uniq % self.n], axis=1)
        ekey = np.array([u * self.n + v for u in range(self.n) for v in nbr[u] if u < v], dtype=np.int64)
        self.pair_is_edge = np.isin(uniq, ekey)

        # link communities (Ahn, Bagrow & Lehmann 2010): two edges sharing a node are
        # linked when the closed-neighbourhood Jaccard of their other endpoints -- the
        # third pair of the wedge they form -- passes the threshold. Per triple, the
        # links between its pair columns (ab-ac via bc, ab-bc via ac, ac-bc via ab),
        # present when the two pairs are edges meeting at the shared node.
        a_, b_, c_ = t[:, 0], t[:, 1], t[:, 2]
        cl = self.closed
        self._links = [((0, 1), cl | (self.centre == a_), self.jpair[:, 2]),
                       ((0, 2), cl | (self.centre == b_), self.jpair[:, 1]),
                       ((1, 2), cl | (self.centre == c_), self.jpair[:, 0])]

    def mask(self, variant: str, t: float) -> np.ndarray:
        if variant in ("hc", "hc_econn"):
            return self.s3 >= t
        if variant == "hc_2test":
            return self.s2 >= t
        if variant == "hc_tri":
            return self.closed & (self.s3 >= t)
        raise ValueError(variant)

    def linkcomm(self, s: float) -> list[set]:
        """Link communities cut at similarity s: components of the graph on edges in
        which two edges sharing a node are joined when the Jaccard of their other two
        endpoints is >= s; each component gives the node set of its edges (kept when
        it has at least three nodes, i.e. at least two edges)."""
        r, c = [], []
        for (x, y), ok, js in self._links:
            sel = ok & (js >= s)
            r.append(self.pair_of[sel, x])
            c.append(self.pair_of[sel, y])
        r, c = np.concatenate(r), np.concatenate(c)
        if len(r) == 0:
            return []
        A = sp.coo_matrix((np.ones(len(r), dtype=np.int8), (r, c)),
                          shape=(self.n_pairs, self.n_pairs)).tocsr()
        _, lab = connected_components(A, directed=False)
        used = np.unique(np.concatenate([r, c]))     # edges with at least one link
        L = np.concatenate([lab[used], lab[used]])
        V = np.concatenate([self.pair_nodes[used, 0], self.pair_nodes[used, 1]])
        key = np.unique(L * self.n + V)
        L, V = key // self.n, key % self.n
        cuts = np.flatnonzero(np.diff(L)) + 1
        return [{self.nodes[i] for i in blk} for blk in np.split(V, cuts) if len(blk) >= 3]

    def linkcomm_density(self, s: float) -> float:
        """Partition density of the link communities cut at s (Ahn, Bagrow & Lehmann
        2010): D = (2/M) sum_c m_c (m_c - n_c + 1) / ((n_c - 2)(n_c - 1)) over edge
        communities with n_c > 2; edges in no link are singletons and contribute 0."""
        r, c = [], []
        for (x, y), ok, js in self._links:
            sel = ok & (js >= s)
            r.append(self.pair_of[sel, x])
            c.append(self.pair_of[sel, y])
        r, c = np.concatenate(r), np.concatenate(c)
        if len(r) == 0 or self.m_edges == 0:
            return 0.0
        A = sp.coo_matrix((np.ones(len(r), dtype=np.int8), (r, c)),
                          shape=(self.n_pairs, self.n_pairs)).tocsr()
        _, lab = connected_components(A, directed=False)
        used = np.unique(np.concatenate([r, c]))
        L = lab[used]
        m_c = np.bincount(L)
        key = np.unique(np.concatenate([L, L]) * self.n
                        + np.concatenate([self.pair_nodes[used, 0], self.pair_nodes[used, 1]]))
        n_c = np.bincount(key // self.n, minlength=len(m_c))
        ok = n_c > 2
        m_c, n_c = m_c[ok].astype(float), n_c[ok].astype(float)
        return float(2.0 / self.m_edges * np.sum(m_c * (m_c - n_c + 1) / ((n_c - 2) * (n_c - 1))))

    def communities(self, variant: str, t: float) -> list[set]:
        keep = np.flatnonzero(self.mask(variant, t))
        if len(keep) == 0:
            return []
        m = len(keep)
        pid = self.pair_of[keep]                     # (m, 3)
        rows = np.repeat(np.arange(m), 3)
        cols = m + pid.ravel()                       # pair nodes after triple nodes
        if variant == "hc_econn":                    # join triples through shared EDGES only
            ok = self.pair_is_edge[pid.ravel()]
            rows, cols = rows[ok], cols[ok]
        A = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                          shape=(m + self.n_pairs, m + self.n_pairs)).tocsr()
        _, lab = connected_components(A, directed=False)
        comp = lab[:m]
        order = np.argsort(comp, kind="stable")
        comp_sorted = comp[order]
        tris = self.triples[keep][order]
        cuts = np.flatnonzero(np.diff(comp_sorted)) + 1
        out = []
        for block in np.split(tris, cuts):
            out.append({self.nodes[i] for i in np.unique(block)})
        return out
