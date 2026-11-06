"""hc_big.BigScores gives the covers of hc_scale.ScaleScores (t_min = 0), the same covers
for every t >= t_min under pruning, and the same overlapping modularity as self_tune."""
import os
import sys

import networkx as nx
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hc_big import BigScores, overlapping_modularity   # noqa: E402
from hc_scale import ScaleScores                       # noqa: E402
from self_tune import overlapping_modularity as om_ref  # noqa: E402

VARIANTS = ["hc", "hc_tri", "hc_2test", "hc_econn"]
TS = [round(0.02 * i, 2) for i in range(1, 31)]


def graphs():
    yield nx.karate_club_graph()
    yield nx.les_miserables_graph()
    for s in range(3):
        yield nx.powerlaw_cluster_graph(400, 4, 0.5, seed=s)
        yield nx.connected_watts_strogatz_graph(300, 8, 0.2, seed=s)
    G = nx.barabasi_albert_graph(600, 3, seed=7)
    G.add_edges_from(nx.star_graph(80).edges())         # a hub with degree-1 leaves
    yield G


def as_edges(G):
    G = nx.convert_node_labels_to_integers(G)
    return G, np.array(G.edges(), dtype=np.int64), G.number_of_nodes()


def sets(node, comm):
    out = {}
    for u, c in zip(node.tolist(), comm.tolist()):
        out.setdefault(c, set()).add(u)
    return sorted(map(frozenset, out.values()), key=sorted)


def canon(cover):
    return sorted(map(frozenset, cover), key=sorted)


def test_equal_scale():
    for G in graphs():
        G, e, n = as_edges(G)
        S = ScaleScores(edges=e, n=n)
        B = BigScores(e, n, block=500, work=2000)
        assert B.n_wedges == S.n_wedges
        for v in VARIANTS:
            for t in TS:
                assert sets(*B.cover(t, v)) == canon(S.communities(t, variant=v)), (v, t)


def test_pruning_exact():
    for G in graphs():
        G, e, n = as_edges(G)
        B0 = BigScores(e, n)
        for tmin in (0.1, 0.2):
            B = BigScores(e, n, t_min=tmin, block=700, work=5000)
            for v in VARIANTS:
                for t in TS:
                    if t >= tmin:
                        assert sets(*B.cover(t, v)) == sets(*B0.cover(t, v)), (tmin, v, t)


def test_modularity():
    for G in graphs():
        G, e, n = as_edges(G)
        A = nx.to_scipy_sparse_array(G, nodelist=range(n), format="csr")
        A = sp.csr_matrix(A)
        idx = {u: u for u in range(n)}
        B = BigScores(e, n)
        for t in (0.1, 0.3):
            node, comm = B.cover(t)
            if len(node) == 0:
                continue
            ref = om_ref(A, idx, sets(node, comm))
            assert abs(overlapping_modularity(A, node, comm) - ref) < 1e-9
