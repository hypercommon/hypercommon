"""
get_node_community must return exactly the global community containing v.

The local expansion and build_hypergraph define the same object: a community is
the union of members over a connected component of the triple-link graph. These
tests pin that equivalence down, across graph families and thresholds.

Regression origin: the local expansion used to track pairs instead of triples.
When node n was admitted from pair (anchor, frontier) it enqueued only
(frontier, n) — never (anchor, n) — so it explored a strict subset of the
admissible triples, and it never checked that the admitted triple had a center.
On ER(50, 0.18, seed=6) at t=0.20 that returned 4 nodes where the global
algorithm returns 45.
"""

import networkx as nx
import pytest

from generators.ring_lattice import ring_lattice
from hypercommon.algorithm import get_communities, get_node_community
from predicates import closed_neighborhood_jaccard_predicate


def assert_equivalent(G, threshold):
    """Every node's local community equals the global component containing it."""
    pred = closed_neighborhood_jaccard_predicate(threshold)
    global_communities = get_communities(G, pred)
    covered = {v for c in global_communities for v in c}

    for v in G.nodes():
        local = get_node_community(G, pred, v)
        candidates = [c for c in global_communities if v in c]

        if v not in covered:
            assert local is None, (
                f"node {v} is in no global community but local returned {local}"
            )
            continue

        assert local is not None, (
            f"node {v} is in a global community of size "
            f"{[len(c) for c in candidates]} but local returned None"
        )
        assert any(local == c for c in candidates), (
            f"node {v} at t={threshold}: local size {len(local)}, "
            f"global candidate sizes {[len(c) for c in candidates]}"
        )


# ================================================================
# The exact case that exposed the bug
# ================================================================

def test_regression_er50_seed6_node1():
    """Local used to return {1,14,29,49}; global returns 45 nodes.

    Node 1 sits in four admissible triples — (1,14,49), (1,29,49), (1,23,29)
    and (1,25,47). The pair-based walk only ever reached the first two, because
    (1,29) was never enqueued and node 23 is the gateway to the rest.
    """
    G = nx.erdos_renyi_graph(50, 0.18, seed=6)
    pred = closed_neighborhood_jaccard_predicate(0.2)

    local = get_node_community(G, pred, 1)
    global_community = next(c for c in get_communities(G, pred) if 1 in c)

    assert local == global_community
    assert len(local) == 45
    # the three nodes the old walk could not reach
    for gateway in (23, 25, 47):
        assert gateway in local


def test_regression_triple_center_not_adjacent_to_both_members():
    """(1,23,29) is admissible with 23 as center — 1 and 29 are not adjacent.

    A purely edge-driven walk cannot see this triple until it steps on 23.
    """
    G = nx.erdos_renyi_graph(50, 0.18, seed=6)
    assert not G.has_edge(1, 29)
    assert G.has_edge(23, 1) and G.has_edge(23, 29)

    local = get_node_community(G, closed_neighborhood_jaccard_predicate(0.2), 1)
    assert {1, 23, 29} <= local


# ================================================================
# Equivalence across graph families
# ================================================================

@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("threshold", [0.1, 0.2, 0.3, 0.4])
def test_equivalence_erdos_renyi(seed, threshold):
    assert_equivalent(nx.erdos_renyi_graph(40, 0.15, seed=seed), threshold)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("threshold", [0.15, 0.25, 0.35])
def test_equivalence_dense_erdos_renyi(seed, threshold):
    assert_equivalent(nx.erdos_renyi_graph(35, 0.35, seed=seed), threshold)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("threshold", [0.1, 0.2, 0.3])
def test_equivalence_barabasi_albert(seed, threshold):
    assert_equivalent(nx.barabasi_albert_graph(45, 3, seed=seed), threshold)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("threshold", [0.1, 0.2, 0.3])
def test_equivalence_watts_strogatz(seed, threshold):
    assert_equivalent(nx.watts_strogatz_graph(45, 6, 0.2, seed=seed), threshold)


@pytest.mark.parametrize("threshold", [0.05, 0.11, 0.2, 0.3])
def test_equivalence_ring_lattice(threshold):
    assert_equivalent(ring_lattice([20] * 3, [6] * 3), threshold)


@pytest.mark.parametrize("threshold", [0.05, 0.11, 0.2, 0.3])
def test_equivalence_uneven_ring_lattice(threshold):
    assert_equivalent(ring_lattice([30, 20, 10], [8, 6, 4]), threshold)


@pytest.mark.parametrize("threshold", [0.05, 0.1, 0.15, 0.2, 0.3, 0.4])
def test_equivalence_karate(threshold):
    assert_equivalent(nx.karate_club_graph(), threshold)


# ================================================================
# Properties of the expansion itself
# ================================================================

def test_result_is_order_independent():
    """Node insertion order must not change the community."""
    G = nx.erdos_renyi_graph(50, 0.18, seed=6)
    pred = closed_neighborhood_jaccard_predicate(0.2)
    expected = get_node_community(G, pred, 1)

    import random

    for k in range(5):
        nodes = list(G.nodes())
        random.Random(k).shuffle(nodes)
        H = nx.Graph()
        H.add_nodes_from(nodes)
        H.add_edges_from(G.edges())
        assert get_node_community(H, pred, 1) == expected


@pytest.mark.parametrize("seed", range(8))
@pytest.mark.parametrize("threshold", [0.15, 0.2, 0.25, 0.3])
def test_a_pair_belongs_to_exactly_one_community(seed, threshold):
    """Nodes overlap between communities; pairs do not.

    Two admissible triples containing the same pair {x, y} share two nodes, so
    build_hypergraph links them and they land in the same component. Hence a
    pair — and therefore an edge — lives in exactly one community. This is what
    makes queued_pairs a sound "process each pair once" rule in expand().

    Note this is a statement about pairs, not about node co-membership: two
    nodes can both appear in two communities without any triple containing both,
    which is a different (and weaker) thing.
    """
    G = nx.erdos_renyi_graph(45, 0.16, seed=seed)
    pred = closed_neighborhood_jaccard_predicate(threshold)

    from hypercommon.hypergraph import build_hypergraph

    H = build_hypergraph(G, pred)
    components = list(nx.connected_components(H))

    pair_to_component = {}
    for index, component in enumerate(components):
        for triple in component:
            a, b, c = triple
            for pair in ((a, b), (a, c), (b, c)):
                previous = pair_to_component.setdefault(pair, index)
                assert previous == index, (
                    f"pair {pair} appears in components {previous} and {index}"
                )


def test_input_graph_is_not_mutated():
    """expand() must not prune edges off the caller's graph."""
    G = nx.erdos_renyi_graph(40, 0.2, seed=3)
    before_edges = set(map(frozenset, G.edges()))
    before_nodes = set(G.nodes())

    get_node_community(G, closed_neighborhood_jaccard_predicate(0.3), 0)

    assert set(map(frozenset, G.edges())) == before_edges
    assert set(G.nodes()) == before_nodes


def test_every_member_is_in_an_admissible_triple():
    """No node may be admitted without a witnessing admissible triple."""
    G = nx.erdos_renyi_graph(45, 0.2, seed=11)
    threshold = 0.2
    pred = closed_neighborhood_jaccard_predicate(threshold)

    from hypercommon.hypernode import HCNode

    hc = {u: HCNode(u, set(G.neighbors(u))) for u in G.nodes()}

    def admissible(a, b, c):
        has_center = (
            (G.has_edge(a, b) and G.has_edge(a, c))
            or (G.has_edge(b, a) and G.has_edge(b, c))
            or (G.has_edge(c, a) and G.has_edge(c, b))
        )
        return has_center and pred(hc[a], hc[b]) and pred(hc[a], hc[c]) and pred(hc[b], hc[c])

    community = get_node_community(G, pred, 0)
    assert community is not None

    for v in community:
        witnessed = any(
            admissible(v, x, y)
            for x in community
            for y in community
            if v != x and v != y and x < y
        )
        assert witnessed, f"node {v} has no admissible triple inside the community"
