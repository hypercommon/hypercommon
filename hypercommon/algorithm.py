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

    # Predicate cache
    predicate_cache = {}

    def check(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in predicate_cache:
            predicate_cache[key] = commonality_predicate(hc[key[0]], hc[key[1]])
        return predicate_cache[key]

    def has_center(a, b, c):
        """True if one of a, b, c is adjacent to the other two."""
        return (
            (G.has_edge(a, b) and G.has_edge(a, c))
            or (G.has_edge(b, a) and G.has_edge(b, c))
            or (G.has_edge(c, a) and G.has_edge(c, b))
        )

    def is_admissible(a, b, c):
        """Same triple condition build_hypergraph uses: a center, plus all three pairs.

        When called from expand(), check(a, b) is redundant — every frontier pair
        came from an already-admitted triple, so its predicate held. It is kept so
        this reads as the full admissibility test rather than something correct
        only under a caller-side invariant; the result is memoised, so it costs a
        dict lookup.
        """
        return has_center(a, b, c) and check(a, b) and check(a, c) and check(b, c)

    def find_initial_triple(v):
        nbrs = list(G.neighbors(v))
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
            for k in G.neighbors(u):
                if k == v:
                    continue
                if check(v, u) and check(u, k) and check(v, k):
                    return v, u, k

        return None

    def expand(initial_triple):
        """Grow the community by walking the same triple-link structure that
        build_hypergraph builds, starting from one admissible triple.

        Frontier entries are unordered pairs {x, y} that belong to some admitted
        triple. A pair is extended by any node n forming another admissible
        triple {x, y, n}; that triple contributes all three of its pairs back to
        the frontier, which is what makes this equivalent to the connected
        component of the hypergraph.
        """
        community = set(initial_triple)
        a, b, c = initial_triple

        def pair_key(x, y):
            return (x, y) if x < y else (y, x)

        potential = deque([pair_key(a, b), pair_key(b, c), pair_key(a, c)])
        queued_pairs = set(potential)

        def push(x, y):
            key = pair_key(x, y)
            if key not in queued_pairs:
                queued_pairs.add(key)
                potential.append(key)

        while potential:
            x, y = potential.popleft()

            # Candidates completing a triple with {x, y} must be adjacent to at
            # least one of them, otherwise no node of the triple can be a center.
            candidates = set(G.neighbors(x))
            candidates |= set(G.neighbors(y))
            candidates.discard(x)
            candidates.discard(y)

            for n in candidates:
                if not is_admissible(x, y, n):
                    continue

                community.add(n)

                # Every pair of the newly admitted triple is a valid frontier —
                # including (x, n) and (y, n).
                push(x, n)
                push(y, n)

        return community

    triple = find_initial_triple(v)
    if triple is None:
        return None

    return expand(triple)