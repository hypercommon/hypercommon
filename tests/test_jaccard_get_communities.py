import networkx as nx
from hypercommon.algorithm import get_communities
from hypercommon.hypernode import HCNode


def test_jaccard_get_communities():
    # --- Build test graph ---
    G = nx.Graph()
    G.add_edges_from([
        (1, 2), (1, 3), (1, 4),
        (2, 3), (2, 4),
        (3, 4),
        (4, 5), (4, 6), (4, 7),
        (5, 6), (5, 7),
        (6, 7)
    ])

    # --- Define Jaccard similarity function ---
    def jaccard(u: HCNode, v: HCNode):
        inter = len(u.neighbors & v.neighbors)
        union = len(u.neighbors | v.neighbors)

        return (inter / union) >= 0.25

    # --- Get communities ---
    communities = get_communities(G, jaccard)

    # --- Expected communities ---
    expected = [
        {1, 2, 3, 4},
        {4, 5, 6, 7}
    ]

    assert communities == expected


test_jaccard_get_communities()