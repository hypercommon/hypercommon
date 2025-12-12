import networkx as nx
from hypercommon.hypergraph import build_hypergraph
from hypercommon.hypernode import HCNode


def test_jaccard_two_components():
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

    # --- Build hypergraph ---
    H = build_hypergraph(G, jaccard)

    # --- Expected hypernodes ---
    expected = {
        (1, 2, 3),
        (1, 2, 4),
        (1, 3, 4),
        (2, 3, 4),
        (4, 5, 6),
        (4, 5, 7),
        (4, 6, 7),
        (5, 6, 7)
    }

    assert set(H.nodes()) == expected

    components = list(nx.connected_components(H))

    assert len(components) == 2


test_jaccard_two_components()