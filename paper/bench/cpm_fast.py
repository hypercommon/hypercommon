"""
Clique percolation (Palla et al. 2005) with the output of networkx's k_clique_communities,
in bounded memory. networkx stores the graph of percolating maximal cliques (pairs sharing
at least k-1 nodes), which on the induced YouTube graph and on the densest small networks
needs more than 100 GB. Here the maximal cliques are enumerated once; the overlaps of a
block of cliques with all others are one sparse product (clique-node incidence times its
transpose), and the percolating pairs of the block are merged into the components found
so far, for every k at once. Only the cliques and one block of overlaps are in memory.

    covers = k_clique_covers(G, ks=(3, 4, 5, 6, 7, 8))   # {k: [frozenset, ...]}
"""
from __future__ import annotations

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

BLOCK_NNZ = 1e8          # bound on the overlaps computed at once (about 1-2 GB)


def k_clique_covers(G, ks=(3, 4, 5, 6, 7, 8), cliques=None):
    ks = sorted(ks)
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    if cliques is None:
        cliques = nx.find_cliques(G)
    rows, cols, sizes = [], [], []
    m = 0
    for c in cliques:
        if len(c) >= ks[0]:
            cols.extend(idx[u] for u in c)
            rows.extend([m] * len(c))
            sizes.append(len(c))
            m += 1
    out = {k: [] for k in ks}
    if m == 0:
        return out
    sizes = np.array(sizes)
    B = sp.csr_matrix((np.ones(len(rows), np.int32), (np.array(rows), np.array(cols))), shape=(m, len(nodes)))
    BT = B.T.tocsr()
    lab = {k: np.arange(m) for k in ks}            # component label of every clique, per k
    # cost of a row: an upper bound on its number of overlapping cliques
    per_node = np.diff(BT.indptr)
    cost = np.add.reduceat(per_node[B.indices], B.indptr[:-1]) if B.nnz else np.zeros(m)
    cum = np.cumsum(cost)
    start = 0
    while start < m:
        stop = int(np.searchsorted(cum, (cum[start - 1] if start else 0) + BLOCK_NNZ, side="right"))
        stop = max(stop, start + 1)
        O = (B[start:stop] @ BT).tocoo()
        r, c, v = O.row + start, O.col, O.data
        keep = c > r
        r, c, v = r[keep], c[keep], v[keep]
        for k in ks:
            sel = (v >= k - 1) & (sizes[r] >= k) & (sizes[c] >= k)
            if not sel.any():
                continue
            L = lab[k]
            # merge: the current components (i -> L[i]) plus the new percolating pairs
            g = sp.coo_matrix((np.ones(m + sel.sum(), np.int8),
                               (np.r_[np.arange(m), r[sel]], np.r_[L, c[sel]])), shape=(m, m))
            _, comp = connected_components(g, directed=False)
            # relabel with the smallest clique index of each component
            first = np.full(comp.max() + 1, m)
            np.minimum.at(first, comp, np.arange(m))
            lab[k] = first[comp]
        start = stop
    for k in ks:
        ok = np.flatnonzero(sizes >= k)
        if not len(ok):
            continue
        Bk = B[ok]
        rr = np.repeat(lab[k][ok], np.diff(Bk.indptr))
        pairs = np.unique(np.c_[rr, Bk.indices], axis=0)
        split = np.flatnonzero(np.diff(pairs[:, 0])) + 1
        out[k] = [frozenset(nodes[j] for j in grp) for grp in np.split(pairs[:, 1], split)]
    return out


if __name__ == "__main__":   # check against networkx on random graphs
    import random
    for s in range(40):
        rnd = random.Random(s)
        G = nx.gnp_random_graph(rnd.randint(20, 150), rnd.uniform(0.05, 0.35), seed=s)
        if s % 2:
            G = nx.relabel_nodes(G, {u: f"v{u}" for u in G})
        got = k_clique_covers(G)
        for k, cov in got.items():
            ref = {frozenset(c) for c in nx.community.k_clique_communities(G, k)}
            assert set(cov) == ref and len(cov) == len(ref), (s, k)
    BLOCK_NNZ = 50      # many small blocks
    for s in range(10):
        G = nx.gnp_random_graph(80, 0.25, seed=s)
        got = k_clique_covers(G)
        for k, cov in got.items():
            assert set(cov) == {frozenset(c) for c in nx.community.k_clique_communities(G, k)}, (s, k)
    print("cpm_fast matches networkx on 50 random graphs")
