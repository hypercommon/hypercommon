"""
Overlap merging for ring-lattice benchmarks.

A ring lattice has disjoint communities, so to benchmark overlapping community
detection we fuse a controlled number of node pairs across rings: for a chosen
pair (u, v) from different rings, v's edges are redirected onto u and v is
removed. The surviving node u then legitimately belongs to both rings, which is
what makes the ground truth overlapping.

Ring sizes may differ, so ring membership comes from the size list rather than
`node // ring_size`.
"""

from __future__ import annotations

import math
import random

import networkx as nx

from .ring_lattice import ring_communities, ring_of_node


def apply_overlap(
        G: nx.Graph,
        sizes: list[int],
        overlap: float,
        rng: random.Random,
) -> dict[int, int]:
    """
    Merge floor(overlap * n) inter-ring node pairs, absorbing v into u.

    Each node takes part in at most one merge, so the requested count is an
    upper bound: the search gives up after k * 20 attempts, which matters when
    overlap is large relative to the number of rings.

    Parameters
    ----------
    G : nx.Graph
        Modified in place — v's edges move to u and v is removed.
    sizes : list[int]
        Ring sizes the graph was built from; defines ring membership.
    overlap : float
        Fraction of n to merge, e.g. 0.05 for 5%.
    rng : random.Random

    Returns
    -------
    dict[int, int]
        {absorbed node: surviving node}
    """
    n = sum(sizes)
    k = math.floor(overlap * n)
    if k <= 0:
        return {}

    node_ring = ring_of_node(sizes)

    merged: dict[int, int] = {}
    used: set[int] = set()

    attempts = 0
    max_attempts = k * 20

    while len(merged) < k and attempts < max_attempts:
        attempts += 1
        u = rng.randrange(n)
        v = rng.randrange(n)

        if u == v or node_ring[u] == node_ring[v]:
            continue
        if u in used or v in used:
            continue
        if not G.has_node(u) or not G.has_node(v):
            continue

        for neighbor in list(G.neighbors(v)):
            if neighbor != u:
                G.add_edge(u, neighbor)
        G.remove_node(v)
        used.add(v)
        used.add(u)
        merged[v] = u

    return merged


def ground_truth_with_overlap(sizes: list[int], merged: dict[int, int]) -> list[set[int]]:
    """
    Ground-truth communities after merging.

    Base case is one community per ring. For each (absorbed v -> surviving u),
    v is dropped from its ring and u takes its place, so u appears in both its
    own ring and v's — the overlap.
    """
    communities = ring_communities(sizes)
    node_ring = ring_of_node(sizes)

    for v, u in merged.items():
        v_ring = node_ring[v]
        communities[v_ring].discard(v)
        communities[v_ring].add(u)

    return communities
