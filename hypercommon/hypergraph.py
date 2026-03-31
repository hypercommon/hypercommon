import networkx as nx
from .hypernode import HCNode


def build_hypergraph(G, commonality_predicate, pair_connection_mode="star"):
    """
    Construct a hypergraph where:

      - Nodes represent unique triplets {i, j, k} where (i,j) and (j,k) are edges in G
        and all three pairs pass the commonality_predicate.
      - Edges connect hypernodes that share exactly two original nodes.
      - commonality_predicate(u, v) is a boolean predicate receiving HCNode objects.

    Returns
    -------
    H : networkx.Graph whose nodes are sorted triplets (i, j, k).
    """

    if not isinstance(G, nx.Graph):
        raise TypeError("G must be a networkx.Graph instance")
    if not callable(commonality_predicate):
        raise TypeError("commonality_predicate must be callable: commonality_predicate(u:HCNode, v:HCNode) -> bool")
    if pair_connection_mode not in ("star", "clique"):
        raise ValueError("pair_connection_mode must be 'star' or 'clique'")

    # ---- Reduce to 2-core (kept for future use) ----
    # G = nx.k_core(G, k=2)

    nodes = list(G.nodes())

    # ---- HCNode cache ----
    hc = {u: HCNode(u, set(G.neighbors(u))) for u in nodes}

    # ---- Hypergraph construction ----
    H = nx.Graph()
    pair_rep = {}
    commonality_cache = {}

    def check(a, b):
        key = (a, b) if a < b else (b, a)
        val = commonality_cache.get(key)
        if val is None:
            val = commonality_predicate(hc[key[0]], hc[key[1]])
            commonality_cache[key] = val
        return val

    def pair_key(x, y):
        return (x, y) if x < y else (y, x)

    seen = set()

    for j in nodes:
        nbrs = list(G.neighbors(j))
        L = len(nbrs)
        for idx_a in range(L):
            i = nbrs[idx_a]
            for idx_b in range(idx_a + 1, L):
                k = nbrs[idx_b]

                triple = tuple(sorted((i, j, k)))
                if triple in seen:
                    continue
                seen.add(triple)

                a, b, c = triple
                if not (check(a, b) and check(a, c) and check(b, c)):
                    continue

                H.add_node(triple, members=triple)

                for x, y in ((a, b), (a, c), (b, c)):
                    key = pair_key(x, y)
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