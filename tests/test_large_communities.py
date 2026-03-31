"""
Tests for two ~10-node communities sharing 2 nodes (9 and 10).

Graph structure:
  Community A (nodes 1-10):
    Core clique: 1-2-3-4
    Periphery: 5 (connects to 1,2,3), 6 (connects to 2,3,4),
               7 (connects to 1,2,3), 8 (connects to 3,4,5)
    Bridge:    9 (connects to 3,4,6), 10 (connects to 2,3,7)

  Community B (nodes 9-18):
    Core clique: 11-12-13-14
    Periphery: 15 (connects to 11,12,13), 16 (connects to 12,13,14),
               17 (connects to 11,12,14), 18 (connects to 13,14,15)
    Bridge:    9  (connects to 11,12,13), 10 (connects to 12,13,14)

At threshold=0.3:
  - All edges within A pass (min J=0.308 for edges 3-9, 3-10)
  - All edges within B pass (min J=0.333 for edges 9-12, 9-13, 10-12, 10-13)
  - No direct edges between A-only and B-only nodes
  - Cross Jaccards (A-only vs B-only) are all < 0.2 -> communities stay separate
  - 170 valid triples in A, 178 valid triples in B
"""

import pytest
import networkx as nx

from hypercommon.algorithm import get_node_community
from predicates import closed_neighborhood_jaccard_predicate


@pytest.fixture
def large_graph():
    G = nx.Graph()
    # Community A
    G.add_edges_from([
        (1, 2), (1, 3),(1, 4), (2, 3), (2, 4), (3, 4),  # core clique
        (1, 5), (2, 5), (3, 5),                     # node 5
        (2, 6), (3, 6), (4, 6),                     # node 6
        (1, 7), (2, 7), (3, 7),                     # node 7
        (3, 8), (4, 8), (5, 8),                     # node 8
        (4, 9), (3, 9), (6, 9),                     # node 9 bridge
        (2, 10), (3, 10), (7, 10),                  # node 10 bridge
    ])
    # Community B
    G.add_edges_from([
        (11,12),(11,13),(11,14),(12,13),(12,14),(13,14),  # core clique
        (11,15),(12,15),(13,15),                           # node 15
        (12,16),(13,16),(14,16),                           # node 16
        (11,17),(12,17),(14,17),                           # node 17
        (13,18),(14,18),(15,18),                           # node 18
        (9,11),(9,12),(9,13),                              # node 9 bridge
        (10,12),(10,13),(10,14),                           # node 10 bridge
    ])
    return G


COMMUNITY_A = set(range(1, 11))
COMMUNITY_B = set(range(9, 19))
A_ONLY = set(range(1, 9))
B_ONLY = set(range(11, 19))
SHARED = {9, 10}
THRESHOLD = 0.3


def test_core_node_A_finds_full_community(large_graph):
    """Node 1 is deep in community A. Should return community containing all A nodes."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 1)
    assert community is not None
    assert A_ONLY <= community, f"A-only nodes missing from community: {A_ONLY - community}"
    assert not (B_ONLY & community), f"B-only nodes leaked into A community: {B_ONLY & community}"


def test_core_node_B_finds_full_community(large_graph):
    """Node 11 is deep in community B. Should return community containing all B nodes."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 11)
    assert community is not None
    assert B_ONLY <= community, f"B-only nodes missing from community: {B_ONLY - community}"
    assert not (A_ONLY & community), f"A-only nodes leaked into B community: {A_ONLY & community}"


def test_periphery_node_A(large_graph):
    """Node 8 is peripheral in A. Should still find community A."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 8)
    assert community is not None
    assert A_ONLY <= community
    assert not (B_ONLY & community)


def test_periphery_node_B(large_graph):
    """Node 18 is peripheral in B. Should still find community B."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 18)
    assert community is not None
    assert B_ONLY <= community
    assert not (A_ONLY & community)


def test_shared_node_9_returns_one_community(large_graph):
    """Node 9 is shared. Must return either community A or B, not a mix."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 9)
    assert community is not None
    in_A = A_ONLY <= community
    in_B = B_ONLY <= community
    assert in_A or in_B, f"Community doesn't match either A or B: {community}"
    assert not (in_A and in_B), "Community contains both A-only and B-only nodes — communities merged"


def test_shared_node_10_returns_one_community(large_graph):
    """Node 10 is shared. Must return either community A or B, not a mix."""
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    community = get_node_community(large_graph, pred, 10)
    assert community is not None
    in_A = A_ONLY <= community
    in_B = B_ONLY <= community
    assert in_A or in_B
    assert not (in_A and in_B)


def test_communities_never_merge(large_graph):
    """
    Cross Jaccards between A-only and B-only nodes are all < 0.2.
    No community should contain both A-only and B-only nodes simultaneously.
    """
    pred = closed_neighborhood_jaccard_predicate(THRESHOLD)
    for v in [1, 2, 3, 11, 12, 13]:
        community = get_node_community(large_graph, pred, v)
        assert community is not None
        has_A_only = bool(A_ONLY & community)
        has_B_only = bool(B_ONLY & community)
        assert not (has_A_only and has_B_only), \
            f"Node {v}: community merged A and B: {community}"


def test_high_threshold_reduces_community(large_graph):
    """
    threshold=0.5: weakest A edges (J~0.308 for 3-9, 3-10) fail.
    Community shrinks but core nodes 1-8 still form valid triples.
    """
    pred = closed_neighborhood_jaccard_predicate(0.5)
    community = get_node_community(large_graph, pred, 1)
    assert community is not None
    # Core clique nodes should still be present
    assert {1, 2, 3, 4} <= community, f"Core nodes missing: {{1,2,3,4}} not in {community}"


def test_no_community_at_very_high_threshold(large_graph):
    """
    threshold=0.9: very few pairs pass. Most nodes return None.
    Only pairs with J>=0.9 are (2,3) J=0.8 and (12,13) J=0.8 — still fail.
    Actually max J in graph is 0.8 so threshold=0.85 -> no valid triples.
    """
    pred = closed_neighborhood_jaccard_predicate(0.85)
    for v in [1, 5, 8, 11, 15, 18]:
        community = get_node_community(large_graph, pred, v)
        assert community is None, f"node {v}: expected None at threshold=0.85, got {community}"
