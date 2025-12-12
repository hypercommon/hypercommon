import networkx as nx
from hypercommon.hypergraph import build_hypergraph
from hypercommon.hypernode import HCNode


def test_two_components_hypergraph():
    # --- Build base graph ---
    G = nx.Graph()
    G.add_edges_from([
        (1, 2), (2, 3), (3, 1),   # component A
        (4, 5), (5, 6), (6, 4)    # component B
    ])

    # --- Define f ---
    def f(u: HCNode, v: HCNode):
        # common neighbors threshold >= 1
        return len(u.neighbors & v.neighbors) >= 1

    # --- Build hypergraph ---
    H = build_hypergraph(G, f)

    # --- Expected hypernodes ---
    expected_hypernodes = {
        (1, 2, 3),
        (4, 5, 6)
    }

    print(H.nodes())

    assert set(H.nodes()) == expected_hypernodes

    # --- Each should be isolated (size=1 component each) ---
    components = list(nx.connected_components(H))
    assert len(components) == 2
    assert {(1, 2, 3)} in components
    assert {(4, 5, 6)} in components


test_two_components_hypergraph()