"""
Edge cases of get_communities / get_node_community found by a stress test against a
brute-force reference (every set of three distinct nodes):

  - a self-loop made a "triple" with a repeated node, so one edge plus a loop gave a
    two-node community;
  - node labels of mixed types (1 and "a", or tuples and ints) crashed the pair order;
  - on a DiGraph the local query and the global cover disagreed.

And the seed triple: at a node in two communities the community returned without a
seed depends on neighbour order; with a seed triple it is the seed's community.
"""
import random

import networkx as nx
import pytest

from hypercommon.algorithm import get_communities, get_node_community
from predicates import closed_neighborhood_jaccard_predicate


def as_sets(cs):
    return sorted(sorted(map(repr, c)) for c in cs)


def two_k5_sharing_node():
    """Two K5 on {0..4} and {4..8}: node 4 is in both communities."""
    G = nx.complete_graph(5)
    G.add_edges_from((a, b) for a in range(4, 9) for b in range(a + 1, 9))
    return G


# --- self-loops are ignored ---

def test_self_loop_with_one_edge_is_no_community():
    G = nx.Graph([(1, 2), (1, 1)])
    pred = closed_neighborhood_jaccard_predicate(0.1)
    assert get_communities(G, pred) == []
    assert get_node_community(G, pred, 1) is None


@pytest.mark.parametrize("t", [0.1, 0.3, 0.5])
def test_self_loops_do_not_change_the_cover(t):
    G = nx.erdos_renyi_graph(40, 0.2, seed=2)
    L = G.copy()
    L.add_edges_from((v, v) for v in range(0, 40, 3))
    pred = closed_neighborhood_jaccard_predicate(t)
    assert as_sets(get_communities(L, pred)) == as_sets(get_communities(G, pred))
    for v in G:
        assert get_node_community(L, pred, v) == get_node_community(G, pred, v)


# --- any hashable labels ---

@pytest.mark.parametrize("t", [0.1, 0.3])
def test_mixed_type_labels(t):
    G = nx.erdos_renyi_graph(30, 0.25, seed=5)
    names = {i: (i if i % 3 == 0 else str(i) if i % 3 == 1 else (i, "x")) for i in G}
    H = nx.relabel_nodes(G, names)
    pred = closed_neighborhood_jaccard_predicate(t)
    want = sorted(sorted(repr(names[v]) for v in c) for c in get_communities(G, pred))
    assert as_sets(get_communities(H, pred)) == want
    glob = [frozenset(c) for c in get_communities(H, pred)]
    for v in H:
        loc = get_node_community(H, pred, v)
        cand = [c for c in glob if v in c]
        assert (loc is None) == (not cand)
        assert loc is None or frozenset(loc) in cand


# --- directed graphs are refused by both ---

def test_directed_graph_rejected():
    D = nx.DiGraph([(0, 1), (1, 2), (2, 0)])
    pred = closed_neighborhood_jaccard_predicate(0.1)
    with pytest.raises(TypeError):
        get_communities(D, pred)
    with pytest.raises(TypeError):
        get_node_community(D, pred, 0)


def test_multigraph_same_as_simple():
    M = nx.MultiGraph([(0, 1), (0, 1), (1, 2), (2, 0), (2, 3), (3, 0)])
    pred = closed_neighborhood_jaccard_predicate(0.1)
    S = nx.Graph(M)
    assert as_sets(get_communities(M, pred)) == as_sets(get_communities(S, pred))
    assert get_node_community(M, pred, 0) == get_node_community(S, pred, 0)


# --- seed triple ---

def test_overlap_node_without_seed_is_order_dependent_but_valid():
    G0 = two_k5_sharing_node()
    pred = closed_neighborhood_jaccard_predicate(0.3)
    comms = [frozenset(c) for c in get_communities(G0, pred)]
    assert sorted(map(sorted, comms)) == [[0, 1, 2, 3, 4], [4, 5, 6, 7, 8]]
    seen = set()
    rs = random.Random(0)
    for _ in range(20):
        G = nx.Graph()
        nodes = list(G0)
        rs.shuffle(nodes)
        G.add_nodes_from(nodes)
        E = list(G0.edges())
        rs.shuffle(E)
        G.add_edges_from(E)
        c = frozenset(get_node_community(G, pred, 4))
        assert c in comms
        seen.add(c)
    assert len(seen) == 2      # both answers occur: hence the seed triple


def test_seed_triple_picks_its_community():
    G = two_k5_sharing_node()
    pred = closed_neighborhood_jaccard_predicate(0.3)
    assert get_node_community(G, pred, 4, triple=(0, 1, 4)) == {0, 1, 2, 3, 4}
    assert get_node_community(G, pred, 4, triple=(4, 7, 8)) == {4, 5, 6, 7, 8}


@pytest.mark.parametrize("seed", range(4))
def test_every_admissible_seed_gives_its_component(seed):
    from hypercommon.hypergraph import build_hypergraph
    G = nx.erdos_renyi_graph(35, 0.2, seed=seed)
    pred = closed_neighborhood_jaccard_predicate(0.2)
    H = build_hypergraph(G, pred)
    for comp in nx.connected_components(H):
        nodes = set().union(*comp)
        for T in comp:
            for v in T:
                assert get_node_community(G, pred, v, triple=T) == nodes


def test_bad_seed_triples():
    G = two_k5_sharing_node()
    pred = closed_neighborhood_jaccard_predicate(0.3)
    with pytest.raises(ValueError):
        get_node_community(G, pred, 4, triple=(0, 1, 2))      # v not in it
    with pytest.raises(ValueError):
        get_node_community(G, pred, 4, triple=(4, 4, 1))      # not distinct
    with pytest.raises(ValueError):
        get_node_community(G, pred, 4, triple=(0, 4, 8))      # pair (0, 8) fails the threshold
    with pytest.raises(nx.NetworkXError):
        get_node_community(G, pred, 4, triple=(4, 1, 99))     # not in G
