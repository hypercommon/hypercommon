"""
HyperCommon for large graphs: the same definition as hc_fast.TripleScores
(test_hc_scale.py checks equality), vectorised with numpy/scipy so that graphs
with tens of millions of wedges fit in memory.

Pipeline
  1. CSR adjacency A (int32 ids), degrees d.
  2. Wedge enumeration per centre j: all neighbour pairs (i, k), i < k, generated
     in vectorised blocks. A closed triangle is reached from all three centres; it
     is kept only from its smallest-id centre among those that see it (i.e. only
     when j < i or the pair (i,k) is not an edge... see _dedup), so each triple is
     scored once.
  3. Pair commonality from C = A @ A (common open neighbours): for a pair (u, v)
     with c = C[u, v],
         adjacent      J = (c + 2) / (d_u + d_v - c)
         non-adjacent  J = c / (d_u + d_v + 2 - c)
     looked up with vectorised sorted-key search, never a Python dict.
  4. Triple score s = min of the three J. For a threshold t, components of the
     bipartite triple--pair incidence graph give the communities.

Memory is O(W) int32/float32 for W wedges; time is dominated by A @ A and by the
wedge blocks, both linear in W up to the cost of the sparse product.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components


def _csr(edges_u, edges_v, n):
    u = np.concatenate([edges_u, edges_v]).astype(np.int64)
    v = np.concatenate([edges_v, edges_u]).astype(np.int64)
    A = sp.csr_matrix((np.ones(len(u), np.int32), (u, v)), shape=(n, n))
    A.sum_duplicates()
    A.data[:] = 1
    A.setdiag(0)
    A.eliminate_zeros()
    A.sort_indices()
    return A


class ScaleScores:
    def __init__(self, G=None, edges=None, n=None, block=2_000_000):
        if G is not None:
            nodes = list(G.nodes())
            idx = {u: i for i, u in enumerate(nodes)}
            e = np.array([(idx[a], idx[b]) for a, b in G.edges()], dtype=np.int64).reshape(-1, 2)
            n = len(nodes)
        else:
            nodes = list(range(n))
            e = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
        self.nodes = nodes
        self.n = n
        A = _csr(e[:, 0], e[:, 1], n)
        self.A = A
        d = np.diff(A.indptr).astype(np.int64)
        self.deg = d
        self.n_wedges = int((d * (d - 1) // 2).sum())

        # edge key set for adjacency tests: key = u * n + v, u < v, sorted
        r = np.repeat(np.arange(n), d)
        c = A.indices.astype(np.int64)
        m = r < c
        self._ekeys = np.sort(r[m] * n + c[m])

        # common-neighbour counts on all distance-<=2 pairs
        C = (A @ A).tocsr()
        C.setdiag(0)
        C.eliminate_zeros()
        Cr = np.repeat(np.arange(n), np.diff(C.indptr))
        Cc = C.indices.astype(np.int64)
        mm = Cr < Cc
        order = np.argsort(Cr[mm] * n + Cc[mm])
        self._ckeys = (Cr[mm] * n + Cc[mm])[order]
        self._cvals = C.data[mm][order].astype(np.int64)

        tri, s3, s2, cl, jp, ce = [], [], [], [], [], []
        for i, j, k, closed in self._wedge_blocks(block):
            a, b, cc = np.sort(np.stack([i, j, k]), axis=0)
            jab, jac, jbc = self._J(a, b), self._J(a, cc), self._J(b, cc)
            sc = np.minimum(np.minimum(jab, jac), jbc)
            # two-test score: an open wedge tests only its two edges, the pairs through the centre j
            inf = np.inf
            e2 = np.minimum(np.minimum(np.where((a == j) | (b == j), jab, inf),
                                       np.where((a == j) | (cc == j), jac, inf)),
                            np.where((b == j) | (cc == j), jbc, inf))
            tri.append(np.stack([a, b, cc], axis=1).astype(np.int32))
            s3.append(sc)
            s2.append(np.where(closed, sc, e2))
            cl.append(closed)
            jp.append(np.stack([jab, jac, jbc], axis=1))
            ce.append(np.where(closed, -1, j).astype(np.int32))
        self.triples = np.concatenate(tri) if tri else np.zeros((0, 3), np.int32)
        self.s3 = np.concatenate(s3) if s3 else np.zeros(0)
        self.s2 = np.concatenate(s2) if s2 else np.zeros(0)
        self.closed = np.concatenate(cl) if cl else np.zeros(0, bool)
        self.jpair = np.concatenate(jp) if jp else np.zeros((0, 3))
        self.centre = np.concatenate(ce) if ce else np.zeros(0, np.int32)

        t = self.triples.astype(np.int64)
        pk = np.concatenate([t[:, 0] * n + t[:, 1], t[:, 0] * n + t[:, 2], t[:, 1] * n + t[:, 2]])
        uniq, inv = np.unique(pk, return_inverse=True)
        self.n_pairs = len(uniq)
        self.pair_of = inv.reshape(3, -1).T.astype(np.int64)
        self.pair_nodes = np.stack([uniq // n, uniq % n], axis=1)
        self.pair_is_edge = self._is_edge(self.pair_nodes[:, 0], self.pair_nodes[:, 1])
        # link communities: links between the pair columns of each triple (see hc_fast)
        cl = self.closed
        self._links = [((0, 1), cl | (self.centre == t[:, 0]), self.jpair[:, 2]),
                       ((0, 2), cl | (self.centre == t[:, 1]), self.jpair[:, 1]),
                       ((1, 2), cl | (self.centre == t[:, 2]), self.jpair[:, 0])]

    # ------------------------------------------------------------------ helpers
    def _is_edge(self, u, v):
        k = np.minimum(u, v) * self.n + np.maximum(u, v)
        pos = np.searchsorted(self._ekeys, k)
        pos = np.minimum(pos, len(self._ekeys) - 1)
        return self._ekeys[pos] == k

    def _common(self, u, v):
        k = np.minimum(u, v) * self.n + np.maximum(u, v)
        pos = np.searchsorted(self._ckeys, k)
        pos = np.minimum(pos, len(self._ckeys) - 1)
        hit = self._ckeys[pos] == k
        return np.where(hit, self._cvals[pos], 0)

    def _J(self, u, v):
        c = self._common(u, v).astype(np.float64)
        du, dv = self.deg[u].astype(np.float64), self.deg[v].astype(np.float64)
        adj = self._is_edge(u, v)
        return np.where(adj, (c + 2) / (du + dv - c), c / (du + dv + 2 - c))

    def _wedge_blocks(self, block):
        """Yield (i, j, k) arrays: j the centre, i < k its neighbours. A closed
        triangle {a, b, c} is seen from each of its three centres; keep only the
        sighting whose centre is the smallest node of the triangle."""
        A, d = self.A, self.deg
        centres = np.flatnonzero(d >= 2)
        cnt = d[centres] * (d[centres] - 1) // 2
        cum = np.cumsum(cnt)
        start = 0
        while start < len(centres):
            stop = int(np.searchsorted(cum, (cum[start - 1] if start else 0) + block, side="right"))
            stop = max(stop, start + 1)
            I, Jc, K = [], [], []
            for j in centres[start:stop]:
                nb = A.indices[A.indptr[j]:A.indptr[j + 1]].astype(np.int64)
                ii, kk = np.triu_indices(len(nb), k=1)
                I.append(nb[ii]); K.append(nb[kk]); Jc.append(np.full(len(ii), j, np.int64))
            i = np.concatenate(I); k = np.concatenate(K); j = np.concatenate(Jc)
            closed = self._is_edge(i, k)
            keep = ~closed | (j < np.minimum(i, k))
            yield i[keep], j[keep], k[keep], closed[keep]
            start = stop

    # ------------------------------------------------------------------ queries
    def linkcomm(self, s):
        """Link communities cut at similarity s (as hc_fast.TripleScores.linkcomm)."""
        r, c = [], []
        for (x, y), ok, js in self._links:
            sel = ok & (js >= s)
            r.append(self.pair_of[sel, x])
            c.append(self.pair_of[sel, y])
        r, c = np.concatenate(r), np.concatenate(c)
        if len(r) == 0:
            return []
        Adj = sp.coo_matrix((np.ones(len(r), np.int8), (r, c)), shape=(self.n_pairs, self.n_pairs)).tocsr()
        _, lab = connected_components(Adj, directed=False)
        used = np.unique(np.concatenate([r, c]))
        L = np.concatenate([lab[used], lab[used]]).astype(np.int64)
        V = np.concatenate([self.pair_nodes[used, 0], self.pair_nodes[used, 1]])
        key = np.unique(L * self.n + V)
        L, V = key // self.n, key % self.n
        cuts = np.flatnonzero(np.diff(L)) + 1
        return [{self.nodes[i] for i in blk} for blk in np.split(V, cuts) if len(blk) >= 3]

    def communities(self, t, triangles_only=False, variant="hc"):
        """variant: hc, hc_2test (edges of open wedges only), hc_tri (triangles only,
        = 3-clique percolation after deleting edges with J < t), hc_econn (hc triples
        joined through shared edges only)."""
        if triangles_only:
            variant = "hc_tri"
        mask = (self.s2 if variant == "hc_2test" else self.s3) >= t
        if variant == "hc_tri":
            mask &= self.closed
        keep = np.flatnonzero(mask)
        if len(keep) == 0:
            return []
        m = len(keep)
        rows = np.repeat(np.arange(m), 3)
        pid = self.pair_of[keep].ravel()
        if variant == "hc_econn":
            ok = self.pair_is_edge[pid]
            rows, pid = rows[ok], pid[ok]
        cols = m + pid
        N = m + self.n_pairs
        Adj = sp.coo_matrix((np.ones(len(rows), np.int8), (rows, cols)), shape=(N, N)).tocsr()
        _, lab = connected_components(Adj, directed=False)
        comp = lab[:m]
        order = np.argsort(comp, kind="stable")
        tris = self.triples[keep][order]
        cuts = np.flatnonzero(np.diff(comp[order])) + 1
        return [{self.nodes[x] for x in np.unique(b)} for b in np.split(tris, cuts)]
