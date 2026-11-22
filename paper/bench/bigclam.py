"""
BigClam (Yang & Leskovec, WSDM 2013): overlapping communities from a nonnegative
affiliation matrix F (nodes x K) fitted by maximum likelihood under

    P(u ~ v) = 1 - exp(-eps0 - F_u . F_v),   eps0 = -ln(1 - 1/n).

A line-for-line port to Python of the authors' reference code (SNAP, snap-adv/agmfast.cpp,
TAGMFast, as run by examples/bigclam with one thread):
  NeighborComInit   seeds: closed neighbourhoods of locally minimal conductance (nodes of
                    degree < 5 never seed; the neighbours of a seed cannot seed); empty
                    communities get 10 random nodes with random weights
  MLEGradAscent     node-by-node projected gradient ascent in a random order, gradient
                    restricted to the neighbours' communities and clipped to [-10, 10],
                    backtracking line search (alpha 0.05, beta 0.3, at most 10 halvings),
                    stop when a pass of 2n updates (after the first 10000) raises the
                    likelihood by at most 1e-4 of its value, or after 1000 n updates
  GetCmtyVV         node u in community c when F_uc >= sqrt(2m / n^2), groups of >= 3 nodes
The parameter is K, the number of communities.
"""
from __future__ import annotations

import numpy as np
from numba import njit


class BigClam:
    def __init__(self, indptr, indices, nodes=None, seed=0):
        """The graph as CSR arrays over 0..n-1 (no self-loops); `nodes` names them."""
        self.indptr, self.indices = np.asarray(indptr, np.int64), np.asarray(indices, np.int64)
        self.n = n = len(self.indptr) - 1
        self.nodes = list(range(n)) if nodes is None else list(nodes)
        self.m = len(self.indices) // 2
        self.eps0 = -np.log(1.0 - 1.0 / n)
        self.rng = np.random.default_rng(seed)
        src = np.repeat(np.arange(n, dtype=np.int64), np.diff(self.indptr))
        up = self.indices > src
        self.ea, self.eb = src[up], self.indices[up]
        self._phi = None

    @classmethod
    def from_nx(cls, G, seed=0):
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        adj = [sorted(idx[v] for v in G.neighbors(u) if v != u) for u in nodes]
        indptr = np.zeros(len(nodes) + 1, dtype=np.int64)
        indptr[1:] = np.cumsum([len(a) for a in adj])
        indices = np.fromiter((v for a in adj for v in a), np.int64, int(indptr[-1]))
        return cls(indptr, indices, nodes, seed)

    def nbrs(self, u):
        return self.indices[self.indptr[u]:self.indptr[u + 1]]

    # -- initialisation -------------------------------------------------------------
    def _conductance_order(self):
        if self._phi is None:
            self._phi = np.argsort(_conductances(self.indptr, self.indices), kind="stable")
        return self._phi

    def init(self, K):
        F = np.zeros((self.n, K))
        invalid = np.zeros(self.n, bool)
        c = 0
        for u in self._conductance_order():
            if c >= K:
                break
            if invalid[u]:
                continue
            F[u, c] = 1.0
            F[self.nbrs(u), c] = 1.0
            invalid[self.nbrs(u)] = True
            c += 1
        for cc in range(K):
            if F[:, cc].sum() == 0.0:
                F[self.rng.integers(0, self.n, 10), cc] = self.rng.random(10)
        self.F, self.sumF = F, F.sum(0)

    # -- fitting (the loop runs in numba: _fit) --------------------------------------------
    def fit(self, thres=1e-4, max_iter_per_node=1000, seed=0):
        _fit(self.indptr, self.indices, self.F, self.sumF, self.eps0, self.ea, self.eb,
             thres, max_iter_per_node * self.n, seed)
        return self

    def communities(self, min_size=3):
        delta = np.sqrt(2.0 * self.m / self.n / self.n)
        out = []
        for c in np.argsort(-self.sumF, kind="stable"):
            if self.sumF[c] < delta:
                continue
            members = np.nonzero(self.F[:, c] >= delta)[0]
            if len(members) >= min_size:
                out.append({self.nodes[i] for i in members})
        return out


def bigclam_covers(G, Ks, seed=0):
    """Yield (K, cover) for every K (K above the number of nodes is skipped)."""
    model = BigClam.from_nx(G, seed)
    for K in Ks:
        if K > model.n:
            continue
        model.rng = np.random.default_rng(seed * 1009 + K)
        model.init(K)
        model.fit(seed=(seed * 1009 + K) % (2**31))
        yield K, model.communities()


@njit(cache=True)
def _conductances(indptr, indices):
    """Conductance of every node's closed neighbourhood (1 for degree < 5), as
    TAGMUtil::GetConductance: cut / min(vol, 2m - vol)."""
    n = len(indptr) - 1
    vol_all = float(len(indices))
    phi = np.ones(n)
    mark = np.full(n, -1, np.int64)
    for u in range(n):
        s, e = indptr[u], indptr[u + 1]
        if e - s < 5:
            continue
        mark[u] = u
        for j in range(s, e):
            mark[indices[j]] = u
        vol = float(e - s)
        cut = 0.0
        for j in range(s, e + 1):
            x = u if j == e else indices[j]
            if j < e:
                vol += indptr[x + 1] - indptr[x]
            for k in range(indptr[x], indptr[x + 1]):
                if mark[indices[k]] != u:
                    cut += 1.0
        den = min(vol, vol_all - vol)
        if den > 0:
            phi[u] = cut / den
        else:
            phi[u] = 0.0 if cut == 0 else 1.0
    return phi


@njit(cache=True)
def _row_lik(indptr, indices, F, sumF, eps0, u, fu):
    K = F.shape[1]
    L = 0.0
    for j in range(indptr[u], indptr[u + 1]):
        v = indices[j]
        dp = 0.0
        for c in range(K):
            dp += F[v, c] * fu[c]
        L += np.log(1.0 - np.exp(-(eps0 + dp))) + dp
    for c in range(K):
        L -= (sumF[c] - F[u, c]) * fu[c]
    return L


@njit(cache=True)
def _likelihood(F, sumF, eps0, ea, eb):
    n, K = F.shape
    L = 0.0
    for i in range(len(ea)):
        dp = 0.0
        for c in range(K):
            dp += F[ea[i], c] * F[eb[i], c]
        L += 2.0 * (np.log(1.0 - np.exp(-(eps0 + dp))) + dp)
    for u in range(n):
        for c in range(K):
            L -= F[u, c] * (sumF[c] - F[u, c])
    return L


@njit(cache=True)
def _fit(indptr, indices, F, sumF, eps0, ea, eb, thres, max_iter, seed):
    """MLEGradAscent of agmfast.cpp (see the module docstring)."""
    np.random.seed(seed)
    n, K = F.shape
    g = np.zeros(K)
    new = np.zeros(K)
    cand = np.zeros(K, np.bool_)
    keep = np.zeros(K, np.bool_)
    it, prev_it = 0, 0
    prev_L, have_prev = 0.0, False
    while it < max_iter:
        order = np.random.permutation(n)
        for u in order:
            it += 1
            s, e = indptr[u], indptr[u + 1]
            if e == s:
                continue
            for c in range(K):
                cand[c] = False
            for j in range(s, e):
                v = indices[j]
                for c in range(K):
                    if F[v, c] > 0.0:
                        cand[c] = True
            any_c = False
            for c in range(K):
                if F[u, c] > 0.0 and not cand[c]:
                    sumF[c] -= F[u, c]
                    F[u, c] = 0.0
                if cand[c]:
                    any_c = True
            if not any_c:
                continue
            for c in range(K):
                g[c] = 0.0
            for j in range(s, e):
                v = indices[j]
                dp = 0.0
                for c in range(K):
                    dp += F[v, c] * F[u, c]
                w = 1.0 / (1.0 - np.exp(-(eps0 + dp)))
                for c in range(K):
                    if cand[c]:
                        g[c] += F[v, c] * w
            gg = 0.0
            for c in range(K):
                keep[c] = False
                if cand[c]:
                    g[c] -= sumF[c] - F[u, c]
                    if not (F[u, c] == 0.0 and g[c] < 0.0) and abs(g[c]) >= 1e-4:
                        keep[c] = True
                        if g[c] > 10.0:
                            g[c] = 10.0
                        elif g[c] < -10.0:
                            g[c] = -10.0
                if not keep[c]:
                    g[c] = 0.0
                gg += g[c] * g[c]
            if gg < 1e-4:
                continue
            L0 = _row_lik(indptr, indices, F, sumF, eps0, u, F[u])
            step = 1.0
            for k in range(10):
                for c in range(K):
                    x = F[u, c] + step * g[c]
                    new[c] = min(max(x, 0.0), 1000.0) if keep[c] else F[u, c]
                if _row_lik(indptr, indices, F, sumF, eps0, u, new) < L0 + 0.05 * step * gg:
                    step *= 0.3
                    if k == 9:
                        step = 0.0
                else:
                    break
            if step == 0.0:
                continue
            for c in range(K):
                if keep[c]:
                    x = F[u, c] + step * g[c]
                    if x <= 0.0:
                        x = 0.0
                    sumF[c] += x - F[u, c]
                    F[u, c] = x
        if it - prev_it >= 2 * n and it > 10000:
            prev_it = it
            L = _likelihood(F, sumF, eps0, ea, eb)
            if have_prev and L - prev_L <= thres * abs(prev_L):
                break
            prev_L, have_prev = L, True
