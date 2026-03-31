"""
Test for overlapping communities sharing exactly one node.

Graph: two triangles sharing node 3.
  Triangle A: 1-2-3 (edges 1-2, 2-3, 1-3)
  Triangle B: 3-4-5 (edges 3-4, 4-5, 3-5)

Jaccard scores at threshold=0.5:
  Triangle A:
    J(1,2) = 1.0   N1={1,2,3},       N2={1,2,3}       inter=3, union=3
    J(1,3) = 0.6   N1={1,2,3},       N3={1,2,3,4,5}   inter=3, union=5
    J(2,3) = 0.6   N2={1,2,3},       N3={1,2,3,4,5}   inter=3, union=5
  Triangle B:
    J(3,4) = 0.6   N3={1,2,3,4,5},   N4={3,4,5}       inter=3, union=5
    J(3,5) = 0.6   N3={1,2,3,4,5},   N5={3,4,5}       inter=3, union=5
    J(4,5) = 1.0   N4={3,4,5},       N5={3,4,5}       inter=3, union=3
  Cross (A vs B, excluding node 3):
    J(1,4) = J(1,5) = J(2,4) = J(2,5) = 0.2  -> all fail at threshold=0.5

Conclusion:
  - Community A = {1,2,3}
  - Community B = {3,4,5}
  - Node 3 is in both (overlapping)
  - The two communities cannot merge since cross pairs fail
"""

import pytest
import networkx as nx

from hypercommon.algorithm import get_node_community
from predicates import closed_neighborhood_jaccard_predicate


@pytest.fixture
def overlapping_graph():
    G = nx.Graph()
    G.add_edges_from([
        (1, 2), (2, 3), (1, 3),  # triangle A
        (3, 4), (4, 5), (3, 5),  # triangle B
    ])
    return G


def test_node1_in_community_A(overlapping_graph):
    """Node 1 is only in triangle A. Community must be {1,2,3}."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(overlapping_graph, pred, 1)
    assert community == {1, 2, 3}


def test_node4_in_community_B(overlapping_graph):
    """Node 4 is only in triangle B. Community must be {3,4,5}."""
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(overlapping_graph, pred, 4)
    assert community == {3, 4, 5}


def test_shared_node3_returns_one_community(overlapping_graph):
    """
    Node 3 is in both communities. get_node_community returns one of them.
    Result must be either {1,2,3} or {3,4,5}.
    """
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(overlapping_graph, pred, 3)
    assert community in [{1, 2, 3}, {3, 4, 5}], \
        f"Expected {{1,2,3}} or {{3,4,5}}, got {community}"


def test_communities_do_not_merge(overlapping_graph):
    """
    Cross pairs J=0.2 < 0.5, so communities cannot merge.
    Neither community should contain nodes from both triangles simultaneously.
    """
    pred = closed_neighborhood_jaccard_predicate(0.5)
    comm1 = get_node_community(overlapping_graph, pred, 1)
    comm4 = get_node_community(overlapping_graph, pred, 4)
    assert not ({1, 2} <= comm4), "Community of node 4 must not contain A-only nodes"
    assert not ({4, 5} <= comm1), "Community of node 1 must not contain B-only nodes"


def test_threshold_too_high_no_community(overlapping_graph):
    """
    threshold=0.7: J(1,3)=J(2,3)=J(3,4)=J(3,5)=0.6 all fail.
    No valid triple exists for any node -> all return None.
    """
    pred = closed_neighborhood_jaccard_predicate(0.7)
    for v in [1, 2, 3, 4, 5]:
        community = get_node_community(overlapping_graph, pred, v)
        assert community is None, f"node {v}: expected None at threshold=0.7, got {community}"
