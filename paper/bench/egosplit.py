"""
Ego-splitting (Epasto, Lattanzi & Paes Leme, KDD 2017): each node is split into personas,
one per cluster of its ego-net (the graph on its neighbours, the node itself removed);
an edge u-v joins u's persona that holds v with v's persona that holds u; a non-overlapping
method clusters the persona graph, and each persona cluster, mapped back to the nodes,
is an overlapping community. A port of the authors' reference code (google-research,
graph_embedding/persona/persona.py: CreatePersonaGraph, PersonaOverlappingClustering)
with its options: local clustering of the ego-nets by connected components or by label
propagation (networkx label_propagation_communities, the reference default); persona
clusters of fewer than 5 personas are dropped (the reference min_cluster_size).

The global clustering of the persona graph is the parameter. The reference offers label
propagation (its default) and greedy modularity; we sweep Leiden (Reichardt-Bornholdt
modularity) over resolutions instead of the single greedy-modularity cut, as for Leiden
itself, and keep label propagation.
"""
from __future__ import annotations

import networkx as nx
from networkx.algorithms.community import label_propagation_communities

MIN_CLUSTER = 5
LOCAL = {"cc": nx.connected_components, "lp": label_propagation_communities}


def persona_graph(G, local="lp"):
    """The persona graph (nodes 0..P-1) and the map persona -> original node."""
    fn = LOCAL[local]
    persona_of = {}                      # (u, v) -> the persona of u that holds neighbour v
    owner = []
    for u in G.nodes():
        nb = [v for v in G.neighbors(u) if v != u]
        ego = G.subgraph(nb)
        for part in fn(ego):
            pid = len(owner)
            owner.append(u)
            for v in part:
                persona_of[(u, v)] = pid
    P = nx.Graph()
    P.add_nodes_from(range(len(owner)))
    for u, v in G.edges():
        if u != v:
            P.add_edge(persona_of[(u, v)], persona_of[(v, u)])
    return P, owner


def _map_back(clusters, owner):
    out = set()
    for c in clusters:
        if len(c) < MIN_CLUSTER:
            continue
        out.add(frozenset(owner[p] for p in c))
    return [set(c) for c in out]


def egosplit_covers(G, grid, seed=0, restarts=3):
    """Yield (param, cover) for each (local, global) of the grid; global is "lp" or a
    Leiden-RB resolution (best of `restarts` runs: 3 as methods2._leiden, 1 on the SNAP graphs
    as big_job.py's Leiden)."""
    import random

    import igraph as ig
    import leidenalg as la
    built = {}
    for local, glob in grid:
        if local not in built:
            P, owner = persona_graph(G, local)
            g = ig.Graph(n=P.number_of_nodes(), edges=list(P.edges()))
            built[local] = (P, owner, g)
        P, owner, g = built[local]
        if glob == "lp":
            random.seed(seed)
            clusters = label_propagation_communities(P)
        else:
            best = None
            for r in range(restarts):
                part = la.find_partition(g, la.RBConfigurationVertexPartition,
                                         resolution_parameter=float(glob), n_iterations=-1,
                                         seed=(seed * 3 + r) % (2**31))
                if best is None or part.quality() > best.quality():
                    best = part
            clusters = [set(c) for c in best if len(c) > 0]
        yield f"{local}:{glob}", _map_back(clusters, owner)
