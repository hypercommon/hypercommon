"""
HyperCommon with a threshold per branch (exploratory; the "first try" of the Discussion).

The covers of HyperCommon at all thresholds t form one hierarchy: a component of the
triple--pair incidence graph at t (triples with score s3 >= t) lies inside one
component at every lower t, and components are disjoint in pairs. Adding the triples
in decreasing order of s3 (Kruskal) gives the single-linkage tree of the pairs, with
lambda = t as the level at which a component splits.

The tree is condensed and clusters are selected as in HDBSCAN (Campello, Moulavit &
Sander 2013; McInnes & Healy 2017), with the pairs as the points:
  * a split at level t is a true split when both sides hold at least `min_size` pairs;
    otherwise the smaller side's pairs fall out of the cluster at t;
  * a cluster born at t_b has stability  sum over its nodes of (t_out - t_b), t_out the
    level at which the node leaves it: when its last pair falls out, or when the cluster
    splits if the node is in a child (unit="pair" counts pairs instead; inside a dense
    ring most pairs leave early, so it favours large merged clusters);
  * excess of mass: going up from the leaves, a cluster is kept when its stability is
    at least the sum of the kept clusters below it; a root that splits is not kept
    (HDBSCAN's allow_single_cluster=False).
A kept cluster's community is the node set of all its pairs at birth, so communities
still overlap in nodes (they are disjoint in pairs). `mode="leaf"` keeps the leaves of
the condensed tree instead.

    TreeHC(ts).communities(min_size, mode="eom")   with ts a hc_fast.TripleScores
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import minimum_spanning_tree


class TreeHC:
    def __init__(self, ts, variant: str = "hc"):
        self.ts = ts
        n = ts.n_pairs
        if variant == "hc":
            keep = np.arange(len(ts.s3))
        elif variant == "hc_tri":
            keep = np.flatnonzero(ts.closed)
        else:
            raise ValueError(variant)
        w = ts.s3[keep]
        po = ts.pair_of[keep]
        # each triple joins its three pairs: edges (p0,p1), (p0,p2) at level w; between two
        # pairs keep the largest level (the edge of the maximum spanning tree)
        a = np.concatenate([po[:, 0], po[:, 0]])
        b = np.concatenate([po[:, 1], po[:, 2]])
        ww = np.concatenate([w, w])
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        key = lo.astype(np.int64) * n + hi
        order = np.lexsort((-ww, key))
        key, ww = key[order], ww[order]
        first = np.r_[True, key[1:] != key[:-1]]
        key, ww = key[first], ww[first]
        # minimum spanning tree on 2 - w  ==  maximum spanning tree on w
        M = sp.coo_matrix((2.0 - ww, (key // n, key % n)), shape=(n, n)).tocsr()
        T = minimum_spanning_tree(M).tocoo()
        lev = 2.0 - T.data
        o = np.argsort(-lev, kind="stable")
        self.n = n
        self.edges = (T.row[o], T.col[o], lev[o])
        # a pair takes part from the level of its best triple
        self.top = np.zeros(n)
        np.maximum.at(self.top, po.ravel(), np.repeat(w, 3))
        self.used = np.zeros(n, bool)
        self.used[po.ravel()] = True
        self._linkage()

    def _linkage(self):
        """Single-linkage merges (child1, child2, level, size) in decreasing level;
        leaves are 0..n-1, merge i is node n+i."""
        n = self.n
        parent = np.arange(n)
        node = np.arange(n)                # union-find root -> current tree node
        size = np.ones(n, dtype=np.int64)

        def find(x):
            r = x
            while parent[r] != r:
                r = parent[r]
            while parent[x] != r:
                parent[x], x = r, parent[x]
            return r

        L = []
        for u, v, lam in zip(*self.edges):
            ru, rv = find(u), find(v)
            if ru == rv:
                continue
            L.append((node[ru], node[rv], lam, size[ru] + size[rv]))
            if size[ru] < size[rv]:
                ru, rv = rv, ru
            parent[rv] = ru
            size[ru] += size[rv]
            node[ru] = n + len(L) - 1
        self.Z = L

    def _condense(self, min_size):
        n, Z = self.n, self.Z
        nz = len(Z)
        left = np.array([z[0] for z in Z], dtype=np.int64)
        right = np.array([z[1] for z in Z], dtype=np.int64)
        lam = np.array([z[2] for z in Z])
        sz = np.ones(n + nz, dtype=np.int64)
        sz[n:] = [z[3] for z in Z]
        has_parent = np.zeros(n + nz, bool)
        has_parent[left] = True
        has_parent[right] = True
        roots = [r for r in range(n, n + nz) if not has_parent[r]]

        def leaves(x):
            out, st = [], [x]
            while st:
                y = st.pop()
                if y < n:
                    out.append(y)
                else:
                    st += [left[y - n], right[y - n]]
            return out

        # condensed tree rows: (parent cluster, child cluster or -1, level, pairs)
        rows_par, rows_child, rows_lam, rows_size, rows_pts = [], [], [], [], []
        birth, cl_root = [], []
        for r in roots:
            if sz[r] < min_size:
                continue
            c = len(birth)
            birth.append(0.0)
            cl_root.append(True)
            st = [(r, c)]
            while st:
                x, c = st.pop()
                if x < n:                       # a single pair left in a cluster
                    rows_par.append(c); rows_child.append(-1)
                    rows_lam.append(self.top[x]); rows_size.append(1); rows_pts.append([x])
                    continue
                a, b, l = left[x - n], right[x - n], lam[x - n]
                big_a, big_b = sz[a] >= min_size, sz[b] >= min_size
                if big_a and big_b:
                    for y in (a, b):
                        k = len(birth)
                        birth.append(l)
                        cl_root.append(False)
                        rows_par.append(c); rows_child.append(k)
                        rows_lam.append(l); rows_size.append(sz[y]); rows_pts.append(None)
                        st.append((y, k))
                else:
                    for y, big in ((a, big_a), (b, big_b)):
                        if big:
                            st.append((y, c))
                        else:
                            pts = leaves(y)
                            rows_par.append(c); rows_child.append(-1)
                            rows_lam.append(l); rows_size.append(len(pts)); rows_pts.append(pts)
        return (np.array(rows_par), np.array(rows_child), np.array(rows_lam),
                np.array(rows_size), rows_pts, np.array(birth), np.array(cl_root))

    def communities(self, min_size: int, mode: str = "eom", unit: str = "node") -> list[set]:
        par, child, lam, size, pts, birth, is_root = self._condense(min_size)
        nc = len(birth)
        if nc == 0:
            return []
        kids = [[] for _ in range(nc)]
        for p, c in zip(par, child):
            if c >= 0:
                kids[p].append(c)
        own = [[] for _ in range(nc)]          # pairs falling out of the cluster itself
        own_lam = [[] for _ in range(nc)]
        for p, c, q, l in zip(par, child, pts, lam):
            if c < 0:
                own[p].extend(q)
                own_lam[p].extend([l] * len(q))
        pn = self.ts.pair_nodes
        # pairs of every cluster at birth (own + descendants), children first
        allp = [None] * nc
        for c in range(nc - 1, -1, -1):
            a = [np.array(own[c], dtype=np.int64)] + [allp[k] for k in kids[c]]
            allp[c] = np.concatenate(a)
        if unit == "pair":
            stab = np.zeros(nc)
            np.add.at(stab, par, (lam - birth[par]) * size)
        else:
            # a node leaves a cluster when its last pair does: when the cluster splits if it
            # is in a child, else at the last level at which one of its pairs falls out
            stab = np.zeros(nc)
            death = np.full(nc, -1.0)
            for p, c, l in zip(par, child, lam):
                if c >= 0:
                    death[p] = l
            for c in range(nc):
                if len(own[c]):
                    q = np.array(own[c], dtype=np.int64)
                    l = np.array(own_lam[c])
                    v = pn[q].ravel()
                    ll = np.repeat(l, 2)
                    ext = {}
                    o = np.lexsort((ll, v))
                    v, ll = v[o], ll[o]
                    last = np.r_[v[1:] != v[:-1], True]
                    ext = dict(zip(v[last].tolist(), ll[last].tolist()))
                else:
                    ext = {}
                if kids[c]:
                    for k in kids[c]:
                        for x in np.unique(pn[allp[k]]).tolist():
                            ext[x] = death[c]
                stab[c] = sum(e - birth[c] for e in ext.values())
        sel = np.zeros(nc, bool)
        # clusters are numbered parents first, so reverse order visits children first
        if mode == "leaf":
            sel = np.array([len(k) == 0 for k in kids])
        else:
            score = stab.copy()
            for c in range(nc - 1, -1, -1):
                if not kids[c]:
                    sel[c] = True
                    continue
                below = sum(score[k] for k in kids[c])
                if stab[c] >= below and not is_root[c]:
                    sel[c] = True
                    score[c] = stab[c]
                else:
                    score[c] = below
            for c in range(nc):                # keep the top-most selected cluster on a path
                if sel[c]:
                    st = list(kids[c])
                    while st:
                        k = st.pop()
                        sel[k] = False
                        st += kids[k]
        out = []
        nodes = self.ts.nodes
        for c in np.flatnonzero(sel):
            v = np.unique(pn[allp[c]])
            if len(v) >= 3:
                out.append({nodes[i] for i in v})
        return out
