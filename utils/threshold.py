from hypercommon.hypernode import HCNode

import networkx as nx


def avg_threshold(
        G,
        commonality_value,
        threshold_multiplier
):
    s = 0
    count = 0

    shortest_path_gen = nx.all_pairs_shortest_path(G, 2)

    hc = {}

    for node in G.nodes():
        hc[node] = HCNode(node, set(G.neighbors(node)))

    for shortest_path in shortest_path_gen:
        node_from = shortest_path[0]
        nodes_to = shortest_path[1]

        del nodes_to[node_from]

        for node_to in nodes_to:
            if node_from < node_to:
                commonality = commonality_value(hc[node_from], hc[node_to])

                s += commonality
                count += 1

    return (s / count) * threshold_multiplier


def representative_thresholds(
        G: nx.Graph,
        commonality_value,
        threshold_multiplier=0.6
):
    u = next(iter(G.nodes()))
    hc_u = HCNode(u, set(G.neighbors(u)))

    neigh1 = set(G.neighbors(u))
    neigh2 = set()

    for v in neigh1:
        neigh2 |= set(G.neighbors(v))
    neigh2.discard(u)

    candidates = neigh1 | neigh2

    values = set()
    for v in candidates:
        hc_v = HCNode(v, set(G.neighbors(v)))
        values.add(commonality_value(hc_u, hc_v))

    values.add(0.0)
    values.add(1.0)

    values.add(avg_threshold(G, commonality_value, threshold_multiplier))
    values = {avg_threshold(G, commonality_value, threshold_multiplier)}

    return sorted(values)