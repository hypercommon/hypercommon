"""
Tests for the ring_lattice generator and its layout helpers.

Rings may differ in both size and degree, so the invariants worth pinning are:
every ring is a separate connected component, every node in ring r has degree
zs[r], and the node ids of ring r are exactly the ground-truth community r.
"""

import networkx as nx
import pytest

from generators.ring_lattice import (
    ring_communities,
    ring_lattice,
    ring_lattice_edge_count,
    ring_of_node,
    ring_offsets,
)


# ================================================================
# Shape and degree
# ================================================================

@pytest.mark.parametrize(
    "sizes,zs",
    [
        ([100] * 4, [16] * 4),
        ([50, 100, 250], [8, 8, 8]),
        ([200, 120, 80], [16, 12, 8]),
        ([30], [4]),
        ([10, 10], [2, 8]),
    ],
)
def test_node_and_edge_counts(sizes, zs):
    G = ring_lattice(sizes, zs)
    assert G.number_of_nodes() == sum(sizes)
    assert G.number_of_edges() == ring_lattice_edge_count(sizes, zs)


@pytest.mark.parametrize(
    "sizes,zs",
    [
        ([100] * 4, [16] * 4),
        ([200, 120, 80], [16, 12, 8]),
        ([40, 25, 15], [10, 6, 4]),
    ],
)
def test_every_node_has_its_rings_degree(sizes, zs):
    G = ring_lattice(sizes, zs)
    for community, z in zip(ring_communities(sizes), zs):
        degrees = {G.degree(v) for v in community}
        assert degrees == {z}


@pytest.mark.parametrize(
    "sizes,zs",
    [
        ([100] * 4, [16] * 4),
        ([200, 120, 80], [16, 12, 8]),
        ([30, 20, 10], [8, 6, 4]),
    ],
)
def test_rings_are_disjoint_components(sizes, zs):
    G = ring_lattice(sizes, zs)
    assert nx.number_connected_components(G) == len(sizes)
    found = {frozenset(c) for c in nx.connected_components(G)}
    expected = {frozenset(c) for c in ring_communities(sizes)}
    assert found == expected


def test_zero_degree_ring_is_isolated_nodes():
    G = ring_lattice([5, 6], [0, 4])
    assert G.number_of_edges() == ring_lattice_edge_count([5, 6], [0, 4]) == 12
    assert all(G.degree(v) == 0 for v in range(5))
    assert all(G.degree(v) == 4 for v in range(5, 11))


# ================================================================
# Layout helpers
# ================================================================

def test_offsets_communities_and_node_mapping_agree():
    sizes = [200, 120, 80]
    offsets = ring_offsets(sizes)
    assert offsets == [0, 200, 320]

    communities = ring_communities(sizes)
    assert [len(c) for c in communities] == sizes
    assert set().union(*communities) == set(range(sum(sizes)))

    mapping = ring_of_node(sizes)
    assert len(mapping) == sum(sizes)
    for r, community in enumerate(communities):
        assert all(mapping[v] == r for v in community)


# ================================================================
# Validation
# ================================================================

def test_mismatched_lengths_rejected():
    with pytest.raises(ValueError, match="same length"):
        ring_lattice([100, 100], [16])


def test_empty_sizes_rejected():
    with pytest.raises(ValueError, match="at least one ring"):
        ring_lattice([], [])


def test_odd_degree_rejected():
    with pytest.raises(ValueError, match="z must be even"):
        ring_lattice([100, 100], [16, 7])


def test_degree_not_smaller_than_ring_rejected():
    """A ring of 30 cannot carry z=32 even when a larger ring in the same
    graph can — the check is per ring, not global."""
    with pytest.raises(ValueError, match="less than the ring size"):
        ring_lattice([200, 30], [16, 32])


def test_nonpositive_size_rejected():
    with pytest.raises(ValueError, match="size must be"):
        ring_lattice([100, 0], [16, 4])


# ================================================================
# Uniform rings still behave like the classic construction
# ================================================================

@pytest.mark.parametrize("n,z,rings", [(400, 16, 4), (500, 8, 10), (60, 6, 3)])
def test_uniform_case_matches_classic_layout(n, z, rings):
    """Equal sizes reproduce the original ring_lattice(n, z, rings) graph."""
    per_ring = n // rings
    G = ring_lattice([per_ring] * rings, [z] * rings)

    expected = nx.Graph()
    expected.add_nodes_from(range(n))
    for r in range(rings):
        offset = r * per_ring
        for i in range(per_ring):
            for k in range(1, z // 2 + 1):
                expected.add_edge(offset + i, offset + (i + k) % per_ring)

    assert set(G.nodes()) == set(expected.nodes())
    assert {frozenset(e) for e in G.edges()} == {frozenset(e) for e in expected.edges()}
