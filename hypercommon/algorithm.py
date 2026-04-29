from collections import deque

import networkx as nx
from hypercommon.hypergraph import build_hypergraph
from hypercommon.hypernode import HCNode

def get_communities(
    G: nx.Graph,
    commonality_predicate,
):
    """
    Compute communities using the Hypercommon method.

    Parameters
    ----------
    G : nx.Graph
        Input graph.
    commonality_predicate : callable
        Function f(u: HCNode, v: HCNode) -> bool.

    Returns
    -------
    list[set]
        List of communities (sets of nodes).
    """

    H = build_hypergraph(G, commonality_predicate)

    communities = []

    for component in nx.connected_components(H):
        nodes = set()
        for hypernode in component:
            nodes.update(hypernode)
        communities.append(nodes)

    return communities


def get_node_community(G: nx.Graph, commonality_predicate, v):
    """
    Find one community containing node v using local hypergraph expansion.

    A community is a set of nodes where every member joined via a valid triple
    (i, j, k) where at least one node is a neighbor of the other two, and all
    three pairs pass the commonality_predicate.

    Parameters
    ----------
    G : nx.Graph
    commonality_predicate : callable (u: HCNode, v: HCNode) -> bool
    v : node in G

    Returns
    -------
    set of nodes, or None if v is not in any community
    """

    # HCNode cache — built from original G, never mutated
    hc = {u: HCNode(u, set(G.neighbors(u))) for u in G.nodes()}

    # Working copy for edge pruning (predicate always uses original hc)
    G_copy = G.copy()

    # Predicate cache
    predicate_cache = {}

    def check(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in predicate_cache:
            predicate_cache[key] = commonality_predicate(hc[key[0]], hc[key[1]])
        return predicate_cache[key]

    def find_initial_triple(v):
        nbrs = list(G_copy.neighbors(v))
        L = len(nbrs)

        # v as center: check pairs among neighbors(v)
        for idx_a in range(L):
            u = nbrs[idx_a]
            for idx_b in range(idx_a + 1, L):
                k = nbrs[idx_b]
                if check(u, v) and check(v, k) and check(u, k):
                    return u, v, k

        # u as center: v -> u -> k
        for u in nbrs:
            for k in G_copy.neighbors(u):
                if k == v:
                    continue
                if check(v, u) and check(u, k) and check(v, k):
                    return v, u, k

        return None

    def expand(initial_triple):
        community = set(initial_triple)
        a, b, c = initial_triple

        # seed potential with all ordered pairs from initial triple
        potential = deque([
            (a, b), (b, a),
            (b, c), (c, b),
            (a, c), (c, a),
        ])

        checked_pairs = set()

        while potential:
            anchor, frontier = potential.popleft()

            if (anchor, frontier) in checked_pairs:
                continue
            checked_pairs.add((anchor, frontier))

            for n in list(G_copy.neighbors(frontier)):
                if (n, frontier) not in checked_pairs:
                    if not check(frontier, n):
                        G_copy.remove_edge(frontier, n)
                        continue

                    if check(anchor, n):
                        if n not in community:
                            community.add(n)
                        potential.append((frontier, n))
                        potential.append((n, frontier))

        return community

    triple = find_initial_triple(v)
    if triple is None:
        return None

    return expand(triple)