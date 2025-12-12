class HCNode:
    """
    A lightweight wrapper around a graph node that exposes:
      - neighbors: set of adjacent nodes
      - degree: number of neighbors

    This is what the user-defined f(u, v) receives.
    """

    __slots__ = ("id", "neighbors", "degree")

    def __init__(self, node_id, neighbors):
        self.id = node_id
        self.neighbors = neighbors
        self.degree = len(neighbors)

    def __repr__(self):
        return f"HCNode({self.id})"