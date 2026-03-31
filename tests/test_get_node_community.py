"""
Tests for get_node_community using closed_neighborhood_jaccard_predicate.

All Jaccard scores are pre-calculated manually and documented per test.
"""

import pytest
import networkx as nx

from hypercommon.algorithm import get_node_community
from predicates import closed_neighborhood_jaccard_predicate


# ================================================================
# GRAPH 1: Simple triangle 1-2-3
# All pairs J=1.0 (closed neighborhoods are identical)
# Any threshold <= 1.0 -> one community {1,2,3}
# ================================================================

@pytest.fixture
def triangle():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (1, 3)])
    return G


def test_triangle_all_nodes_one_community(triangle):
    """J(1,2)=J(1,3)=J(2,3)=1.0. threshold=0.5 -> community={1,2,3}"""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    for v in [1, 2, 3]:
        community = get_node_community(triangle, pred, v)
        assert community == {1, 2, 3}, f"node {v}: expected {{1,2,3}}, got {community}"


def test_triangle_threshold_too_high(triangle):
    """
    Triangle J=1.0 for all pairs, so even threshold=1.0 should work.
    """
    pred = closed_neighborhood_jaccard_predicate(1.0)
    community = get_node_community(triangle, pred, 1)
    assert community == {1, 2, 3}


def test_triangle_node_not_in_graph(triangle):
    """Node 99 not in graph — should raise."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    with pytest.raises(Exception):
        get_node_community(triangle, pred, 99)


# ================================================================
# GRAPH 2: Two triangles sharing edge 2-3
# edges: 1-2, 2-3, 1-3, 3-4, 2-4
#
# Jaccard scores:
#   J(1,2) = 0.75  (N1={1,2,3}, N2={1,2,3,4}, inter=3, union=4)
#   J(1,3) = 0.75  (N1={1,2,3}, N3={1,2,3,4}, inter=3, union=4)
#   J(2,3) = 1.0   (N2=N3={1,2,3,4})
#   J(2,4) = 0.75  (N2={1,2,3,4}, N4={2,3,4}, inter=3, union=4)
#   J(3,4) = 0.75  (N3={1,2,3,4}, N4={2,3,4}, inter=3, union=4)
#   J(1,4) = 0.5   (N1={1,2,3}, N4={2,3,4}, inter=2, union=4)
# ================================================================

@pytest.fixture
def two_triangles():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (1, 3), (3, 4), (2, 4)])
    return G


def test_two_triangles_low_threshold_one_community(two_triangles):
    """
    threshold=0.5: all pairs >= 0.5, so all four nodes in one community.
    J(1,4)=0.5 is the weakest link — still passes.
    """
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(two_triangles, pred, 1)
    assert community == {1, 2, 3, 4}


def test_two_triangles_mid_threshold_node1_community(two_triangles):
    """
    threshold=0.6: J(1,4)=0.5 fails.
    Node 1 can still join via triple (1,2,3) since J(1,2)=J(1,3)=0.75, J(2,3)=1.0.
    Node 4 can join via triple (2,3,4) since J(2,4)=J(3,4)=0.75, J(2,3)=1.0.
    But can 1 and 4 be in same community? Yes — connected through 2 and 3.
    Community = {1,2,3,4}.
    """
    pred = closed_neighborhood_jaccard_predicate(0.6)
    community = get_node_community(two_triangles, pred, 1)
    assert community == {1, 2, 3, 4}


def test_two_triangles_high_threshold_node1_excluded(two_triangles):
    """
    threshold=0.8: J(1,2)=J(1,3)=0.75 fail.
    Node 1 cannot form any valid triple -> returns None.
    Node 4 also cannot (J(2,4)=J(3,4)=0.75 fail).
    Only nodes 2 and 3 pass all checks but they need a third -> None.
    """
    pred = closed_neighborhood_jaccard_predicate(0.8)
    community = get_node_community(two_triangles, pred, 1)
    assert community is None


def test_two_triangles_high_threshold_center_nodes(two_triangles):
    """
    threshold=0.8: J(2,3)=1.0 passes but no valid triple exists
    since all other pairs fail. Node 2 returns None.
    """
    pred = closed_neighborhood_jaccard_predicate(0.8)
    community = get_node_community(two_triangles, pred, 2)
    assert community is None


# ================================================================
# GRAPH 3: Two separate triangles
# Triangle A: 1-2-3, Triangle B: 4-5-6, no edges between them
#
# Jaccard scores within each triangle = 1.0 (same as graph 1)
# Jaccard cross triangles = 0.0 (no shared neighbors)
# ================================================================

@pytest.fixture
def two_separate_triangles():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (1, 3)])  # triangle A
    G.add_edges_from([(4, 5), (5, 6), (4, 6)])  # triangle B
    return G


def test_separate_triangles_node_in_A(two_separate_triangles):
    """Node 1 in triangle A -> community {1,2,3}, never reaches triangle B."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(two_separate_triangles, pred, 1)
    assert community == {1, 2, 3}


def test_separate_triangles_node_in_B(two_separate_triangles):
    """Node 4 in triangle B -> community {4,5,6}, never reaches triangle A."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(two_separate_triangles, pred, 4)
    assert community == {4, 5, 6}


def test_separate_triangles_communities_are_disjoint(two_separate_triangles):
    """Communities of node 1 and node 4 must be disjoint."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    comm_A = get_node_community(two_separate_triangles, pred, 1)
    comm_B = get_node_community(two_separate_triangles, pred, 4)
    assert comm_A == {1, 2, 3}
    assert comm_B == {4, 5, 6}
    assert comm_A.isdisjoint(comm_B)


# ================================================================
# GRAPH 4: K4 clique — all 4 nodes fully connected
# All pairs J=1.0 (identical closed neighborhoods)
# All thresholds <= 1.0 -> one community {1,2,3,4}
# ================================================================

@pytest.fixture
def k4():
    G = nx.Graph()
    G.add_edges_from([
        (1, 2), (1, 3), (1, 4),
        (2, 3), (2, 4),
        (3, 4)
    ])
    return G


def test_k4_full_community(k4):
    """All pairs J=1.0 -> entire clique is one community."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    for v in [1, 2, 3, 4]:
        community = get_node_community(k4, pred, v)
        assert community == {1, 2, 3, 4}, f"node {v}: expected {{1,2,3,4}}, got {community}"


def test_k4_threshold_1(k4):
    """threshold=1.0 still works since all J=1.0."""
    pred = closed_neighborhood_jaccard_predicate(1.0)
    community = get_node_community(k4, pred, 1)
    assert community == {1, 2, 3, 4}


# ================================================================
# GRAPH 5: Isolated node
# Node 7 has no edges -> no valid triple -> None
# ================================================================

def test_isolated_node_returns_none():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (1, 3)])
    G.add_node(7)
    pred = closed_neighborhood_jaccard_predicate(0.0)
    community = get_node_community(G, pred, 7)
    assert community is None


# ================================================================
# GRAPH 6: Node reachable only via wedge (not triangle)
# edges: 1-2, 2-3 (path, no edge 1-3)
#
# Jaccard:
#   J(1,2): N1={1,2}, N2={1,2,3}, inter=2, union=3 = 0.667
#   J(2,3): N2={1,2,3}, N3={2,3}, inter=2, union=3 = 0.667
#   J(1,3): N1={1,2}, N3={2,3}, inter=1, union=3 = 0.333
# triple (1,2,3): check(1,2)=0.667, check(2,3)=0.667, check(1,3)=0.333
# ================================================================

def test_wedge_low_threshold():
    """
    threshold=0.3: J(1,3)=0.333 >= 0.3 -> triple (1,2,3) valid -> community {1,2,3}.
    """
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    pred = closed_neighborhood_jaccard_predicate(0.3)
    community = get_node_community(G, pred, 1)
    assert community == {1, 2, 3}


def test_wedge_high_threshold():
    """
    threshold=0.4: J(1,3)=0.333 < 0.4 -> triple (1,2,3) fails -> None.
    """
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    pred = closed_neighborhood_jaccard_predicate(0.4)
    community = get_node_community(G, pred, 1)
    assert community is None
