"""
HyperCommon for very large graphs (SNAP YouTube): the covers of hc_scale.ScaleScores
for the variants hc, hc_tri, hc_2test and hc_econn (test_hc_big.py checks equality),
with three changes that keep memory near O(W) small integers for W scored wedges.

  1. Degree-ratio pruning, exact for every t >= t_min. Since
         J(u, v) <= (min(d_u, d_v) + 1) / (max(d_u, d_v) + 1),
     an edge whose ratio is below t_min fails every test at t >= t_min, and so does
     every triple containing it. Wedges are enumerated on the kept edges only; the
     degrees, the closed/open test and the common-neighbour counts use the full graph.
  2. Common-neighbour counts are computed only for the pairs of the scored triples,
     from row blocks of A @ A, never the whole product.
  3. Components are taken on the pair graph (each admitted triple links its pairs),
     and covers are returned as membership arrays rather than Python sets.

Link communities are not provided: their similarity J(i, k) does not involve the two
edges' own scores, so the edge pruning above does not apply to them.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components


def _csr(u, v, n):
    A = sp.csr_matrix((np.ones(2 * len(u), np.int32), (np.r_[u, v], np.r_[v, u])), shape=(n, n))
    A.sum_duplicates()
    A.data[:] = 1
    A.setdiag(0)
    A.eliminate_zeros()
    A.sort_indices()
    return A


def _member(keys, sorted_keys):
    pos = np.minimum(np.searchsorted(sorted_keys, keys), len(sorted_keys) - 1)
    return sorted_keys[pos] == keys


class BigScores:
    def __init__(self, edges, n, t_min=0.0, block=4_000_000, work=150_000_000, log=None):
        log = log or (lambda *a: None)
        e = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
        self.n = n
        A = _csr(e[:, 0], e[:, 1], n)
        self.A = A
        d = np.diff(A.indptr).astype(np.int64)
        self.deg = d
        self.n_wedges_full = int((d * (d - 1) // 2).sum())

        r = np.repeat(np.arange(n), d)
        c = A.indices.astype(np.int64)
        up = r < c
        eu, ev = r[up], c[up]
        del r, c, up
        self._ekeys = eu * n + ev                        # sorted: CSR rows, sorted columns
        if t_min > 0:
            ok = (np.minimum(d[eu], d[ev]) + 1) / (np.maximum(d[eu], d[ev]) + 1) >= t_min
            self._kkeys = self._ekeys[ok]
            Ak = _csr(eu[ok], ev[ok], n)
        else:
            self._kkeys = self._ekeys
            Ak = A
        del eu, ev
        dk = np.diff(Ak.indptr).astype(np.int64)
        self.n_wedges = int((dk * (dk - 1) // 2).sum())
        log(f"edges {len(self._ekeys)}, kept {len(self._kkeys)}; wedges {self.n_wedges_full} -> {self.n_wedges}")

        # 1. wedges on the kept edges
        tri, cen = [], []
        for i, j, k, closed in self._wedge_blocks(Ak, dk, block):
            t = np.sort(np.stack([i, j, k]), axis=0).T.astype(np.int32)
            tri.append(t)
            cen.append(np.where(closed, -1, j).astype(np.int32))
        del Ak
        self.triples = np.concatenate(tri) if tri else np.zeros((0, 3), np.int32)
        del tri
        self.centre = np.concatenate(cen) if cen else np.zeros(0, np.int32)
        del cen
        self.closed = self.centre < 0
        W = len(self.triples)
        log(f"triples scored {W}")

        # 2. pairs of the triples
        t = self.triples
        pk = np.empty(3 * W, np.int64)
        pk[:W] = t[:, 0].astype(np.int64) * n + t[:, 1]
        pk[W:2 * W] = t[:, 0].astype(np.int64) * n + t[:, 2]
        pk[2 * W:] = t[:, 1].astype(np.int64) * n + t[:, 2]
        uniq, inv = np.unique(pk, return_inverse=True)
        del pk
        self.n_pairs = len(uniq)
        self.pair_of = np.ascontiguousarray(inv.astype(np.int32 if len(uniq) < 2**31 else np.int64).reshape(3, -1).T)
        del inv
        pu, pv = uniq // n, uniq % n
        self.pair_is_edge = _member(uniq, self._ekeys)
        log(f"pairs {self.n_pairs}")

        # 3. common neighbours of those pairs, from row blocks of A @ A
        cnt = self._common(uniq, pu, work)
        du, dv = d[pu].astype(np.float64), d[pv].astype(np.float64)
        cf = cnt.astype(np.float64)
        Jp = np.where(self.pair_is_edge, (cf + 2) / (du + dv - cf), cf / (du + dv + 2 - cf))
        del du, dv, cf, cnt
        self.pair_nodes = np.stack([pu, pv], axis=1).astype(np.int32)
        del uniq, pu, pv

        # 4. triple scores: all three pairs (s3); for open wedges, the two edges only (s2)
        P = self.pair_of
        j0, j1, j2 = Jp[P[:, 0]], Jp[P[:, 1]], Jp[P[:, 2]]
        self.s3 = np.minimum(np.minimum(j0, j1), j2)
        cj = self.centre
        inf = np.inf
        # pairs (a,b), (a,c), (b,c) contain the centre j when j is one of their nodes
        e2 = np.minimum(np.minimum(np.where((t[:, 0] == cj) | (t[:, 1] == cj), j0, inf),
                                   np.where((t[:, 0] == cj) | (t[:, 2] == cj), j1, inf)),
                        np.where((t[:, 1] == cj) | (t[:, 2] == cj), j2, inf))
        self.s2 = np.where(self.closed, self.s3, e2)
        del j0, j1, j2, e2, Jp

    # ------------------------------------------------------------------ helpers
    def _wedge_blocks(self, Ak, dk, block):
        """(i, j, k): centre j, i < k its kept neighbours. Closed triangles are tested
        on the full graph; one with a pruned edge fails every test and is dropped; one
        with three kept edges is kept from its smallest-id centre only."""
        centres = np.flatnonzero(dk >= 2)
        cum = np.cumsum(dk[centres] * (dk[centres] - 1) // 2)
        start = 0
        while start < len(centres):
            stop = int(np.searchsorted(cum, (cum[start - 1] if start else 0) + block, side="right"))
            stop = max(stop, start + 1)
            I, Jc, K = [], [], []
            for j in centres[start:stop]:
                nb = Ak.indices[Ak.indptr[j]:Ak.indptr[j + 1]].astype(np.int64)
                ii, kk = np.triu_indices(len(nb), k=1)
                I.append(nb[ii]); K.append(nb[kk]); Jc.append(np.full(len(ii), j, np.int64))
            i = np.concatenate(I); k = np.concatenate(K); j = np.concatenate(Jc)
            key = np.minimum(i, k) * self.n + np.maximum(i, k)
            closed = _member(key, self._ekeys)
            keep = ~closed | (j < np.minimum(i, k))
            if self._kkeys is not self._ekeys:
                keep &= ~closed | _member(key, self._kkeys)
            yield i[keep], j[keep], k[keep], closed[keep]
            start = stop

    def _common(self, keys, pu, work):
        """|N(u) & N(v)| for the sorted pair keys u * n + v (u < v), from row blocks of
        A @ A sized so that each block costs about `work` multiply-adds."""
        A, n = self.A, self.n
        out = np.zeros(len(keys), np.int64)
        rows = np.unique(pu)
        cost = A @ self.deg.astype(np.float64)              # work of row u in A @ A
        cum = np.cumsum(cost[rows])
        start = 0
        while start < len(rows):
            stop = int(np.searchsorted(cum, (cum[start - 1] if start else 0) + work, side="right"))
            stop = max(stop, start + 1)
            R = rows[start:stop]
            S = (A[R] @ A).tocsr()
            S.sort_indices()
            sk = np.repeat(R.astype(np.int64), np.diff(S.indptr)) * n + S.indices
            lo, hi = np.searchsorted(keys, [int(R[0]) * n, (int(R[-1]) + 1) * n])
            q = keys[lo:hi]
            pos = np.minimum(np.searchsorted(sk, q), max(len(sk) - 1, 0))
            if len(sk):
                hit = sk[pos] == q
                out[lo:hi] = np.where(hit, S.data[pos], 0)
            start = stop
        return out

    # ------------------------------------------------------------------ queries
    def cover(self, t, variant="hc"):
        """Cover at threshold t as membership arrays (node, community), community ids
        0..C-1, sorted by community then node."""
        mask = (self.s2 if variant == "hc_2test" else self.s3) >= t
        if variant == "hc_tri":
            mask &= self.closed
        keep = np.flatnonzero(mask)
        if len(keep) == 0:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        P = self.pair_of[keep]
        if variant == "hc_econn":
            # join through edges only: an open wedge's non-adjacent pair carries no link
            E = self.pair_is_edge[P]
            # first and second edge pair of each triple (every admitted triple has >= 2)
            order = np.argsort(~E, axis=1, kind="stable")
            p0 = np.take_along_axis(P, order[:, :1], 1).ravel()
            p1 = np.take_along_axis(P, order[:, 1:2], 1).ravel()
            p2 = np.where(E.all(axis=1), np.take_along_axis(P, order[:, 2:3], 1).ravel(), p0)
        else:
            p0, p1, p2 = P[:, 0], P[:, 1], P[:, 2]
        r = np.concatenate([p0, p0]).astype(np.int64)
        c = np.concatenate([p1, p2]).astype(np.int64)
        Adj = sp.coo_matrix((np.ones(len(r), np.int8), (r, c)), shape=(self.n_pairs, self.n_pairs)).tocsr()
        del r, c
        _, lab = connected_components(Adj, directed=False)
        del Adj
        comp = lab[p0].astype(np.int64)
        tri = self.triples[keep]
        key = np.unique(np.concatenate([comp, comp, comp]) * self.n
                        + np.concatenate([tri[:, 0], tri[:, 1], tri[:, 2]]).astype(np.int64))
        cc, node = key // self.n, key % self.n
        _, cc = np.unique(cc, return_inverse=True)
        return node, cc.ravel()


def overlapping_modularity(A, node, comm):
    """self_tune.overlapping_modularity for a cover given as membership arrays."""
    n = A.shape[0]
    M2 = A.nnz
    if M2 == 0 or len(node) == 0:
        return 0.0
    k = np.asarray(A.sum(axis=1)).ravel().astype(np.float64)
    O = np.bincount(node, minlength=n).astype(np.float64)
    C = int(comm.max()) + 1
    Wm = sp.csr_matrix((1.0 / O[node], (node, comm)), shape=(n, C))
    Af = A.astype(np.float64)
    q1 = float((Af @ Wm).multiply(Wm).sum())
    q2 = float((np.asarray(Wm.T @ k).ravel() ** 2).sum())
    return (q1 - q2 / M2) / M2
