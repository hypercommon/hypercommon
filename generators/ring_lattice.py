import networkx as nx


def ring_lattice(
        n: int,
        z: int,
        rings: int
) -> nx.Graph:
    """
    Generate a ring lattice (or multiple disjoint ring lattices).

    Parameters
    ----------
    n : int
        Total number of nodes.
    z : int
        Even degree of each node (z/2 neighbors on each side).
    rings : int
        Number of identical disjoint rings.

    Returns
    -------
    nx.Graph
    """

    if rings < 1:
        raise ValueError("rings must be >= 1")

    if n % rings != 0:
        raise ValueError("n must be divisible by rings")

    if z % 2 != 0:
        raise ValueError("z must be even")

    n_per_ring = n // rings

    if z >= n_per_ring:
        raise ValueError("z must be less than nodes per ring")

    G = nx.Graph()
    G.add_nodes_from(range(n))

    half = z // 2

    for r in range(rings):
        offset = r * n_per_ring

        for i in range(n_per_ring):
            u = offset + i
            for k in range(1, half + 1):
                v = offset + (i + k) % n_per_ring
                G.add_edge(u, v)

    return G