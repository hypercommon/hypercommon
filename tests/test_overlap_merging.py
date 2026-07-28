"""
Tests for apply_overlap / ground_truth_with_overlap.

Ring membership is derived from the size list, so these hold for unequal rings
as well as the uniform case the experiments used previously.
"""

import random

import pytest

from generators.overlap import apply_overlap, ground_truth_with_overlap
from generators.ring_lattice import ring_lattice, ring_of_node


def build(sizes, zs):
    return ring_lattice(sizes, zs)


# ================================================================
# Merge counts and structure
# ================================================================

@pytest.mark.parametrize(
    "sizes,zs,overlap,expected",
    [
        ([100] * 4, [16] * 4, 0.0, 0),
        ([100] * 4, [16] * 4, 0.01, 4),
        ([100] * 4, [16] * 4, 0.05, 20),
        ([200, 120, 80], [16, 12, 8], 0.05, 20),
        ([200, 120, 80], [16, 12, 8], 0.10, 40),
    ],
)
def test_merge_count(sizes, zs, overlap, expected):
    G = build(sizes, zs)
    merged = apply_overlap(G, sizes, overlap, random.Random(42))
    assert len(merged) == expected


def test_absorbed_nodes_are_removed_from_the_graph():
    sizes, zs = [200, 120, 80], [16, 12, 8]
    G = build(sizes, zs)
    before = G.number_of_nodes()

    merged = apply_overlap(G, sizes, 0.05, random.Random(1))

    assert G.number_of_nodes() == before - len(merged)
    for absorbed, surviving in merged.items():
        assert not G.has_node(absorbed)
        assert G.has_node(surviving)


def test_merges_only_ever_cross_rings():
    """A merge inside one ring would not create overlap, so it must not happen."""
    sizes, zs = [200, 120, 80], [16, 12, 8]
    G = build(sizes, zs)
    node_ring = ring_of_node(sizes)

    merged = apply_overlap(G, sizes, 0.08, random.Random(3))

    assert merged
    for absorbed, surviving in merged.items():
        assert node_ring[absorbed] != node_ring[surviving]


def test_each_node_participates_in_at_most_one_merge():
    sizes, zs = [150, 150, 100], [12, 12, 8]
    G = build(sizes, zs)

    merged = apply_overlap(G, sizes, 0.10, random.Random(5))

    involved = list(merged.keys()) + list(merged.values())
    assert len(involved) == len(set(involved))


def test_surviving_node_inherits_the_absorbed_neighbourhood():
    sizes, zs = [60, 40], [6, 4]
    G = build(sizes, zs)
    original = {v: set(G.neighbors(v)) for v in G.nodes()}

    merged = apply_overlap(G, sizes, 0.10, random.Random(11))

    for absorbed, surviving in merged.items():
        for neighbour in original[absorbed]:
            if neighbour != surviving and G.has_node(neighbour):
                assert G.has_edge(surviving, neighbour)


# ================================================================
# Ground truth
# ================================================================

def test_ground_truth_places_survivor_in_both_rings():
    sizes, zs = [200, 120, 80], [16, 12, 8]
    G = build(sizes, zs)
    node_ring = ring_of_node(sizes)

    merged = apply_overlap(G, sizes, 0.05, random.Random(7))
    truth = ground_truth_with_overlap(sizes, merged)

    assert len(truth) == len(sizes)
    for absorbed, surviving in merged.items():
        assert absorbed not in truth[node_ring[absorbed]]
        assert surviving in truth[node_ring[absorbed]]
        assert surviving in truth[node_ring[surviving]]


def test_ground_truth_without_overlap_is_the_plain_rings():
    sizes = [200, 120, 80]
    truth = ground_truth_with_overlap(sizes, {})
    assert [len(c) for c in truth] == sizes
    assert set().union(*truth) == set(range(sum(sizes)))


def test_ground_truth_total_membership_is_conserved():
    """Each merge moves one member; it never changes the total count."""
    sizes, zs = [150, 150, 100], [12, 12, 8]
    G = build(sizes, zs)

    merged = apply_overlap(G, sizes, 0.06, random.Random(13))
    truth = ground_truth_with_overlap(sizes, merged)

    assert sum(len(c) for c in truth) == sum(sizes)


# ================================================================
# Determinism
# ================================================================

def test_same_seed_gives_the_same_merge():
    sizes, zs = [200, 120, 80], [16, 12, 8]
    a = apply_overlap(build(sizes, zs), sizes, 0.05, random.Random(99))
    b = apply_overlap(build(sizes, zs), sizes, 0.05, random.Random(99))
    assert a == b


def test_different_seeds_give_different_merges():
    sizes, zs = [200, 120, 80], [16, 12, 8]
    a = apply_overlap(build(sizes, zs), sizes, 0.05, random.Random(1))
    b = apply_overlap(build(sizes, zs), sizes, 0.05, random.Random(2))
    assert a != b
