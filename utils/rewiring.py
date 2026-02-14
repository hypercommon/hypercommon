import random
import networkx as nx


def rewire_one_edge(
        G: nx.Graph,
        keep_node: int,
        replace_node: int,
        rng: random.Random
) -> None:
    # pick a new endpoint w
    candidates = set(G.nodes())
    candidates.discard(keep_node)
    candidates -= set(G.neighbors(keep_node))  # prevents self-loop and multi-edge

    if not candidates:
        print('Error finding candidate edge')
        return  # extremely rare / degenerate

    new_node = rng.choice(tuple(candidates))

    G.remove_edge(keep_node, replace_node)
    G.add_edge(keep_node, new_node)


def rewire_step(
        G: nx.Graph,
        edge_stack: list[tuple[int, int]],
        k: int,
        rng: random.Random
) -> None:
    for _ in range(k):
        u, v = edge_stack.pop()

        if random.getrandbits(1):
            keep_node, replace_node = u, v
        else:
            keep_node, replace_node = v, u

        rewire_one_edge(G, keep_node, replace_node, rng)
