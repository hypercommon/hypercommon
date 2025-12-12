import networkx as nx
from .hypernode import HCNode


def build_hypergraph(G, f):
    """
    Construct a hypergraph where:

      - Nodes represent triplets (i, j, k)
      - Edges connect triplets sharing exactly two original nodes
      - f(u, v) is a boolean function receiving HCNode objects

    Returns
    -------
    H : networkx.Graph
        Graph whose nodes are triplets (i,j,k).
    """

    # ---- Basic validation ----
    if not isinstance(G, nx.Graph):
        raise TypeError("G must be a networkx.Graph instance")
    if not callable(f):
        raise TypeError("f must be callable f(u:HCNode, v:HCNode) -> bool")

    # ---- Reduce to 2-core ----
    G = nx.k_core(G, k=2)
    nodes = sorted(G.nodes())
    if len(nodes) < 3:
        return nx.Graph()

    # ---- Two-step neighborhoods ----
    two_step = {}
    hc = {}

    for u, lengths in nx.all_pairs_shortest_path_length(G, cutoff=2):
        neighbors_u = set()
        bigger = []

        for v, dist in lengths.items():
            if dist == 1:
                neighbors_u.add(v)

            if v > u:
                bigger.append(v)

        bigger.sort()

        hc[u] = HCNode(u, neighbors_u)
        two_step[u] = bigger

    # ---- Hypergraph construction ----
    H = nx.Graph()
    pair_to_hypernodes = {}

    f_cache = {}

    def check(a, b):
        # Invariant: a < b always holds by construction (ordered iteration).
        # Nodes are indexed by integers; (a, b) uniquely identifies the unordered pair.
        key = (a, b)

        if key not in f_cache:
            f_cache[key] = f(hc[a], hc[b])

        return f_cache[key]

    for i in nodes:
        neigh_i = two_step.get(i, [])
        L = len(neigh_i)

        for idx_j in range(L):
            j = neigh_i[idx_j]
            neigh_j = two_step.get(j, [])

            for idx_k in range(idx_j + 1, L):
                k = neigh_i[idx_k]

                # structural pruning
                if k not in neigh_j:
                    continue

                # semantic pruning using f(HCNode, HCNode)
                if not (check(i, j) and check(i, k) and check(j, k)):
                    continue

                # valid hypernode
                hypernode = (i, j, k)
                H.add_node(hypernode, members=hypernode)

                # connect via shared pairs
                for key in [(i, j), (i, k), (j, k)]:
                    if key in pair_to_hypernodes:
                        for other in pair_to_hypernodes[key]:
                            H.add_edge(hypernode, other)
                        pair_to_hypernodes[key].append(hypernode)
                    else:
                        pair_to_hypernodes[key] = [hypernode]

    return H