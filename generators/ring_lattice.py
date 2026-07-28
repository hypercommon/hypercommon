import networkx as nx


def ring_lattice(
        sizes: list[int],
        zs: list[int],
) -> nx.Graph:
    """
    Generate a set of disjoint ring lattices, one per entry of `sizes`.

    Ring r has sizes[r] nodes, each joined to its zs[r]/2 nearest neighbours on
    either side. Rings are laid out consecutively, so ring r owns the node ids
    [offset_r, offset_r + sizes[r]) where offset_r is the sum of the preceding
    sizes. Total node count is sum(sizes).

    Parameters
    ----------
    sizes : list[int]
        Nodes per ring. Length is the number of rings.
    zs : list[int]
        Even degree for each ring; must be the same length as `sizes`.
        A ring of size s requires zs[r] < s.

    Returns
    -------
    nx.Graph

    Examples
    --------
    Four equal rings of 100 at degree 16:

    >>> G = ring_lattice([100] * 4, [16] * 4)

    Unequal rings, each with its own degree:

    >>> G = ring_lattice([200, 120, 80], [16, 12, 8])
    """

    if len(sizes) != len(zs):
        raise ValueError(f"sizes and zs must be the same length: {len(sizes)} != {len(zs)}")

    if not sizes:
        raise ValueError("sizes must contain at least one ring")

    for r, (size, z) in enumerate(zip(sizes, zs)):
        if size < 1:
            raise ValueError(f"ring {r}: size must be >= 1, got {size}")
        if z % 2 != 0:
            raise ValueError(f"ring {r}: z must be even, got {z}")
        if z < 0:
            raise ValueError(f"ring {r}: z must be >= 0, got {z}")
        if z >= size:
            raise ValueError(f"ring {r}: z must be less than the ring size, got z={z} size={size}")

    G = nx.Graph()
    G.add_nodes_from(range(sum(sizes)))

    offset = 0
    for size, z in zip(sizes, zs):
        half = z // 2
        for i in range(size):
            u = offset + i
            for k in range(1, half + 1):
                G.add_edge(u, offset + (i + k) % size)
        offset += size

    return G


def ring_offsets(sizes: list[int]) -> list[int]:
    """First node id of each ring, i.e. the cumulative sum of `sizes`."""
    offsets = []
    running = 0
    for size in sizes:
        offsets.append(running)
        running += size
    return offsets


def ring_of_node(sizes: list[int]) -> list[int]:
    """Map node id -> ring index, for the layout `ring_lattice` produces."""
    mapping = []
    for r, size in enumerate(sizes):
        mapping.extend([r] * size)
    return mapping


def ring_communities(sizes: list[int]) -> list[set[int]]:
    """Ground-truth communities: one set of node ids per ring."""
    communities = []
    offset = 0
    for size in sizes:
        communities.append(set(range(offset, offset + size)))
        offset += size
    return communities


def ring_lattice_edge_count(sizes: list[int], zs: list[int]) -> int:
    """Edges in `ring_lattice(sizes, zs)`, without building the graph."""
    return sum(size * z // 2 for size, z in zip(sizes, zs))
