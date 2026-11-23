import networkx as nx
from .hypernode import HCNode


def _ordered(a, b):
    """Canonical order of a node pair. Nodes of one comparable type keep their natural
    order; mixed types (say 1 and "a") fall back to (type name, repr), so any hashable
    labels work."""
    try:
        return (a, b) if a < b else (b, a)
    except TypeError:
        return (a, b) if (type(a).__name__, repr(a)) < (type(b).__name__, repr(b)) else (b, a)


def _sorted_triple(i, j, k):
    try:
        return tuple(sorted((i, j, k)))
    except TypeError:
        return tuple(sorted((i, j, k), key=lambda x: (type(x).__name__, repr(x))))


def _neighbors(G, u):
    """Neighbours of u without u itself: a self-loop is not a pair of distinct nodes."""
    return [w for w in G.neighbors(u) if w != u]


def _check_graph(G):
    if not isinstance(G, nx.Graph):
        raise TypeError("G must be a networkx.Graph instance")
    if G.is_directed():
        raise TypeError("G must be undirected (use G.to_undirected())")


def build_hypergraph(G, commonality_predicate, pair_connection_mode="star"):
    """
    Construct a hypergraph where:

      - Nodes represent unique triplets {i, j, k} of distinct nodes where at least one
        node is a neighbor of the other two (i.e. there exists a center node adjacent
        to both others), and all three pairs pass the commonality_predicate.
      - Edges connect hypernodes that share exactly two original nodes.
      - commonality_predicate(u, v) is a boolean predicate receiving HCNode objects;
        it should be symmetric. Self-loops are ignored.

    Returns
    -------
    H : networkx.Graph whose nodes are sorted triplets (i, j, k).
    """

    _check_graph(G)
    if not callable(commonality_predicate):
        raise TypeError("commonality_predicate must be callable: commonality_predicate(u:HCNode, v:HCNode) -> bool")
    if pair_connection_mode not in ("star", "clique"):
        raise ValueError("pair_connection_mode must be 'star' or 'clique'")

    # ---- Reduce to 2-core (kept for future use) ----
    # G = nx.k_core(G, k=2)

    nodes = list(G.nodes())

    # ---- HCNode cache ----
    hc = {u: HCNode(u, set(_neighbors(G, u))) for u in nodes}

    # ---- Hypergraph construction ----
    H = nx.Graph()
    pair_rep = {}
    commonality_cache = {}

    def check(a, b):
        key = _ordered(a, b)
        val = commonality_cache.get(key)
        if val is None:
            val = commonality_predicate(hc[key[0]], hc[key[1]])
            commonality_cache[key] = val
        return val

    seen = set()

    for j in nodes:
        nbrs = _neighbors(G, j)
        L = len(nbrs)
        for idx_a in range(L):
            i = nbrs[idx_a]
            for idx_b in range(idx_a + 1, L):
                k = nbrs[idx_b]

                triple = _sorted_triple(i, j, k)
                if triple in seen:
                    continue
                seen.add(triple)

                a, b, c = triple
                if not (check(a, b) and check(a, c) and check(b, c)):
                    continue

                H.add_node(triple, members=triple)

                for x, y in ((a, b), (a, c), (b, c)):
                    key = _ordered(x, y)
                    rep = pair_rep.get(key)
                    if rep is None:
                        if pair_connection_mode == "star":
                            pair_rep[key] = triple
                        else:
                            pair_rep[key] = [triple]
                    else:
                        if pair_connection_mode == "star":
                            H.add_edge(triple, rep)
                        else:
                            for node in pair_rep[key]:
                                H.add_edge(triple, node)
                            pair_rep[key].append(triple)

    return H
