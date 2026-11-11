"""
Every method under one interface: `run_method(name, G, ts, seed)` returns a list
of (param, cover) for that method's whole predeclared parameter grid.

The grids are fixed here, before any result is seen, and every parameterised
method is tuned the same two ways (see run.py): frozen from calibration graphs,
and per-graph oracle. Parameter-free methods have a one-point grid.

Stochastic methods are seeded from the graph seed so reruns are identical.
"""

from __future__ import annotations

import random
import warnings

import numpy as np

warnings.filterwarnings("ignore")

T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]          # 0.02 .. 0.60

GRIDS = {
    "hc": T_GRID,
    "hc_2test": T_GRID,
    "hc_tri": T_GRID,                       # = Jaccard-filtered 3-clique percolation
    "cpm": [3, 4, 5],
    "leiden": [0.25, 0.5, 1.0, 2.0, 4.0, 8.0],   # RB-configuration resolution
    "infomap": [None],
    "lpa": [None],
    "walktrap": [None],
    "slpa": [0.05, 0.1, 0.2, 0.3, 0.45],
    "demon": [0.1, 0.25, 0.5, 0.75],
    "angel": [0.1, 0.25, 0.5, 0.75],
    "scan": [(e, m) for e in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8) for m in (2, 3, 5)],
    "scan_j": [(e, m) for e in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6) for m in (2, 3, 5)],
    "linkcomm": [None],
}

# Grid extensions, added after a boundary check found these methods' best parameter
# at the edge of the original grid on many graphs. Run as separate methods ("<name>_x")
# and merged into their base method by analysis.py.
GRIDS_X = {
    "slpa": [0.01, 0.02],
    "demon": [0.9, 1.0],
    "angel": [0.9, 1.0],
    "cpm": [6, 7],
    "scan": [(e, m) for e in (0.1, 0.2) for m in (2, 3, 5)],
}

HC_VARIANTS = ("hc", "hc_2test", "hc_tri")
SLOW = {"walktrap"}    # skipped when n is large (see run.py)


def _sets(nc):
    return [set(c) for c in nc.communities if len(c) > 0]


def _structural_similarity(G, cosine=True):
    """SCAN similarity on closed neighbourhoods, per edge: cosine (Xu et al. 2007) or Jaccard."""
    nb = {u: set(G.adj[u]) | {u} for u in G}
    sim = {}
    for u, v in G.edges():
        a, b = nb[u], nb[v]
        inter = len(a & b)
        s = inter / (len(a) * len(b)) ** 0.5 if cosine else inter / (len(a) + len(b) - inter)
        sim[(u, v)] = sim[(v, u)] = s
    return sim


def _scan(G, sim, eps, mu):
    """
    SCAN: a core has >= mu members in its eps-neighbourhood (itself included);
    clusters are cores joined by eps-similar edges, plus every non-core node in a
    core's eps-neighbourhood (attached to each such cluster -- a node bordering
    two clusters keeps both, the overlapping reading of SCAN's border nodes).
    Hubs and outliers stay unassigned.
    """
    eps_nb = {u: [v for v in G.adj[u] if sim[(u, v)] >= eps] for u in G}
    core = {u for u in G if len(eps_nb[u]) + 1 >= mu}
    parent = {u: u for u in core}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u in core:
        for v in eps_nb[u]:
            if v in core:
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
    clusters = {}
    for u in core:
        r = find(u)
        clusters.setdefault(r, set()).add(u)
        clusters[r].update(eps_nb[u])
    return [c for c in clusters.values() if len(c) >= 3]


def run_method(name, G, ts, seed):
    """Yield (param, cover). `ts` is the precomputed TripleScores of G."""
    if name == "slpa100":                          # SLPA with the original paper's 100 iterations
        from cdlib import algorithms as A
        random.seed(seed)
        np.random.seed(seed % (2**32))
        for r in GRIDS["slpa"] + GRIDS_X["slpa"]:
            yield r, _sets(A.slpa(G, t=100, r=r))
        return
    if name.endswith("_x"):                        # grid extension of a base method
        base = name[:-2]
        saved = GRIDS[base]
        GRIDS[base] = GRIDS_X[base]
        try:
            yield from run_method(base, G, ts, seed)
        finally:
            GRIDS[base] = saved
        return
    if name.startswith("hcf_"):                    # HyperCommon with commonality f, quantile grid
        from hc_family import Q_GRID
        f = name[4:]
        for q in Q_GRID:
            yield q, ts.communities_f(f, q)
        return
    if name in HC_VARIANTS:
        for t in GRIDS[name]:
            yield t, ts.communities(name, t)
        return

    from cdlib import algorithms as A
    random.seed(seed)
    np.random.seed(seed % (2**32))

    if name == "cpm":
        # the output of cdlib's kclique (networkx k_clique_communities), in bounded memory
        from cpm_fast import k_clique_covers
        covers = k_clique_covers(G, GRIDS[name])
        for k in GRIDS[name]:
            yield k, [set(c) for c in covers[k]]
    elif name == "leiden":
        import igraph as ig
        import leidenalg as la
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        g = ig.Graph(n=len(nodes), edges=[(idx[u], idx[v]) for u, v in G.edges()])
        for gamma in GRIDS[name]:
            part = la.find_partition(g, la.RBConfigurationVertexPartition,
                                     resolution_parameter=gamma, seed=seed % (2**31))
            yield gamma, [{nodes[i] for i in c} for c in part if len(c) > 0]
    elif name == "infomap":
        from infomap import Infomap
        im = Infomap(f"--two-level --silent --seed {seed % (2**31) + 1}")
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        for u, v in G.edges():
            im.add_link(idx[u], idx[v])
        im.run()
        groups = {}
        for node in im.tree:
            if node.is_leaf:
                groups.setdefault(node.module_id, set()).add(nodes[node.node_id])
        yield None, [c for c in groups.values() if len(c) > 0]
    elif name == "lpa":
        yield None, _sets(A.label_propagation(G))
    elif name == "walktrap":
        yield None, _sets(A.walktrap(G))
    elif name == "slpa":
        for r in GRIDS[name]:
            yield r, _sets(A.slpa(G, t=21, r=r))
    elif name == "demon":
        for e in GRIDS[name]:
            yield e, _sets(A.demon(G, epsilon=e, min_com_size=3))
    elif name == "angel":
        for th in GRIDS[name]:
            yield th, _sets(A.angel(G, threshold=th, min_community_size=3))
    elif name in ("scan", "scan_j"):
        sim = _structural_similarity(G, cosine=(name == "scan"))
        for e, m in GRIDS[name]:
            yield (e, m), _scan(G, sim, e, m)
    elif name == "linkcomm":
        ec = A.hierarchical_link_community(G)
        covers = []
        for edges in ec.communities:
            nodes = set()
            for u, v in edges:
                nodes.add(u)
                nodes.add(v)
            if len(nodes) >= 3:
                covers.append(nodes)
        yield None, covers
    else:
        raise ValueError(name)
