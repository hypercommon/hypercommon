import networkx as nx
from hypercommon.hypergraph import build_hypergraph


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