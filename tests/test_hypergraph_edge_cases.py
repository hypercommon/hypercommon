"""
Edge-case tests for the hypergraph algorithm (build_hypergraph).

These tests cover cases that could break the code or expose incorrect behavior.
Run with: python -m pytest tests/test_hypergraph_edge_cases.py -v
Or run with: python -c "from tests.test_hypergraph_edge_cases import *; run_all()"
"""
import networkx as nx
from hypercommon.hypergraph import build_hypergraph
from hypercommon.hypernode import HCNode

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False


def _always_true(_u: HCNode, _v: HCNode) -> bool:
    return True


def _always_false(_u: HCNode, _v: HCNode) -> bool:
    return False


# --- Empty and small graphs (should not crash) ---


def test_empty_graph():
    """Empty graph: no nodes, iterator yields nothing; should return empty H."""
    G = nx.Graph()
    H = build_hypergraph(G, _always_true)
    assert H.number_of_nodes() == 0
    assert H.number_of_edges() == 0


def test_single_node():
    """Single node: two_step[node] = []; no (j,k) pairs; should return empty H."""
    G = nx.Graph()
    G.add_node(0)
    H = build_hypergraph(G, _always_true)
    assert H.number_of_nodes() == 0
    assert H.number_of_edges() == 0


def test_two_nodes_one_edge():
    """Two nodes: at most one node in two_step[i]; no idx_j, idx_k with j < k; empty H."""
    G = nx.Graph()
    G.add_edge(0, 1)
    H = build_hypergraph(G, _always_true)
    assert H.number_of_nodes() == 0
    assert H.number_of_edges() == 0


# --- Triangle vs wedge (structural pruning) ---


def test_wedge_produces_no_hypernode():
    """Wedge 1-2-3 (edges 1-2, 2-3 only): must NOT create hypernode (1,2,3)."""
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    H = build_hypergraph(G, _always_true)
    # Algorithm should restrict to triangles; wedge is not a triangle.
    assert H.number_of_nodes() == 0, "Wedge (1-2-3) must not form a hypernode"
    assert H.number_of_edges() == 0


def test_triangle_produces_hypernode():
    """Triangle 1-2-3: must create exactly one hypernode (1,2,3)."""
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1)])
    H = build_hypergraph(G, _always_true)
    assert set(H.nodes()) == {(1, 2, 3)}
    assert H.number_of_edges() == 0


# --- Predicate raises ---


def test_predicate_raises_propagates():
    """If commonality_predicate raises for a pair, build_hypergraph should raise."""
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1)])

    def bad_predicate(u: HCNode, v: HCNode) -> bool:
        if u.id == 1 and v.id == 2:
            raise ValueError("predicate failed")
        return True

    if HAS_PYTEST:
        with pytest.raises(ValueError, match="predicate failed"):
            build_hypergraph(G, bad_predicate)
    else:
        try:
            build_hypergraph(G, bad_predicate)
            assert False, "expected ValueError"
        except ValueError as e:
            assert "predicate failed" in str(e)


# --- Node types: sorted() requires comparable nodes ---


def test_non_comparable_nodes_break_sorted():
    """Nodes that are not mutually comparable break sorted(G.nodes())."""
    G = nx.Graph()
    G.add_nodes_from([1, "a", (3,)])  # int, str, tuple
    G.add_edges_from([(1, "a"), ("a", (3,)), ((3,), 1)])  # triangle-ish

    if HAS_PYTEST:
        with pytest.raises(TypeError):
            build_hypergraph(G, _always_true)
    else:
        try:
            build_hypergraph(G, _always_true)
            assert False, "expected TypeError from sorted()"
        except TypeError:
            pass


# --- Type validation (already in code) ---


def test_non_graph_raises():
    """Passing non-Graph (e.g. None or list) raises TypeError."""
    # Use something that is definitely not nx.Graph (DiGraph may be accepted in some NetworkX versions)
    if HAS_PYTEST:
        with pytest.raises(TypeError, match="networkx.Graph"):
            build_hypergraph(None, _always_true)
    else:
        try:
            build_hypergraph(None, _always_true)
            assert False, "expected TypeError"
        except TypeError as e:
            assert "networkx.Graph" in str(e) or "Graph" in str(e)


def test_non_callable_predicate_raises():
    """Passing non-callable predicate raises TypeError."""
    G = nx.Graph()
    G.add_edge(1, 2)

    if HAS_PYTEST:
        with pytest.raises(TypeError, match="callable"):
            build_hypergraph(G, "not a function")
    else:
        try:
            build_hypergraph(G, "not a function")
            assert False, "expected TypeError"
        except TypeError as e:
            assert "callable" in str(e).lower()


def run_all():
    """Run all edge-case tests without pytest."""
    test_empty_graph()
    test_single_node()
    test_two_nodes_one_edge()
    test_wedge_produces_no_hypernode()
    test_triangle_produces_hypernode()
    test_predicate_raises_propagates()
    test_non_comparable_nodes_break_sorted()
    test_non_graph_raises()
    test_non_callable_predicate_raises()
    print("All edge-case tests passed.")


if __name__ == "__main__":
    run_all()
