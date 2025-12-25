import networkx as nx
from hypercommon.hypergraph import build_hypergraph


def get_communities(
    G: nx.Graph,
    similarity_fn,
):
    """
    Compute communities using the Hypercommon method.

    Parameters
    ----------
    G : nx.Graph
        Input graph.
    similarity_fn : callable
        Function f(u: HCNode, v: HCNode) -> bool.

    Returns
    -------
    list[set]
        List of communities (sets of nodes).
    """

    H = build_hypergraph(G, similarity_fn)

    communities = []

    for component in nx.connected_components(H):
        nodes = set()
        for hypernode in component:
            nodes.update(hypernode)
        communities.append(nodes)

    return communities