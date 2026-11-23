"""
get_communities and get_node_community against the definition itself, by brute force:
every set of three distinct nodes with a centre whose three pairs pass the predicate is
admitted, admitted triples sharing a pair are joined, a community is the node union of a
class. Small graphs of many kinds, several predicates.
"""
import itertools
from collections import defaultdict

import networkx as nx
import pytest

from generators.ring_lattice import ring_lattice
from hypercommon.algorithm import get_communities, get_node_community
from hypercommon.hypernode import HCNode
from predicates import closed_neighborhood_jaccard_predicate


def brute(G, pred):
    nodes = list(G)
    hc = {u: HCNode(u, set(G.neighbors(u)) - {u}) for u in nodes}
    adj = lambda a, b: a != b and G.has_edge(a, b)
    parent, owner, admitted = {}, {}, []

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, c in itertools.combinations(nodes, 3):
        if not ((adj(a, b) and adj(a, c)) or (adj(b, a) and adj(b, c)) or (adj(c, a) and adj(c, b))):
            continue
        if not all(pred(hc[x], hc[y]) for x, y in ((a, b), (a, c), (b, c))):
            continue
        T = frozenset((a, b, c))
        parent[T] = T
        admitted.append(T)
        for p in map(frozenset, itertools.combinations(T, 2)):
            if p in owner:
                ra, rb = find(T), find(owner[p])
                if ra != rb:
                    parent[ra] = rb
            else:
                owner[p] = T
    comps = defaultdict(set)
    for T in admitted:
        comps[find(T)] |= T
    return [frozenset(c) for c in comps.values()]


def open_jaccard(t):
    def p(u, v):
        U, V = u.neighbors - {u.id}, v.neighbors - {v.id}
        return (len(U & V) / len(U | V) if U | V else 0) >= t
    return p


GRAPHS = {
    "K6": nx.complete_graph(6),
    "star": nx.star_graph(7),
    "wheel": nx.wheel_graph(9),
    "grid (tuple nodes)": nx.grid_2d_graph(4, 5),
    "petersen": nx.petersen_graph(),
    "karate": nx.karate_club_graph(),
    "florentine (str nodes)": nx.florentine_families_graph(),
    "barbell": nx.barbell_graph(5, 2),
    "windmill": nx.windmill_graph(4, 4),
    "caveman": nx.relaxed_caveman_graph(4, 6, 0.15, seed=1),
    "disjoint union": nx.disjoint_union_all([nx.complete_graph(4), nx.cycle_graph(5), nx.empty_graph(2)]),
    "ring lattice": ring_lattice([20, 12], [6, 4]),
    "ER sparse": nx.erdos_renyi_graph(30, 0.15, seed=1),
    "ER dense": nx.erdos_renyi_graph(25, 0.4, seed=2),
    "BA": nx.barabasi_albert_graph(30, 3, seed=3),
    "geometric": nx.random_geometric_graph(30, 0.3, seed=4),
}
PREDS = {
    "J>=0": closed_neighborhood_jaccard_predicate(0.0),
    "J>=.15": closed_neighborhood_jaccard_predicate(0.15),
    "J>=.3": closed_neighborhood_jaccard_predicate(0.3),
    "J>=1": closed_neighborhood_jaccard_predicate(1.0),
    "open J>=.2": open_jaccard(0.2),
    "never": lambda u, v: False,
}


@pytest.mark.parametrize("pname", PREDS)
@pytest.mark.parametrize("gname", GRAPHS)
def test_matches_definition(gname, pname):
    G, pred = GRAPHS[gname], PREDS[pname]
    ref = brute(G, pred)
    glob = [frozenset(c) for c in get_communities(G, pred)]
    key = lambda cs: sorted(sorted(map(repr, c)) for c in cs)
    assert key(glob) == key(ref)
    for v in G:
        loc = get_node_community(G, pred, v)
        cand = [c for c in ref if v in c]
        assert (loc is None) == (not cand)
        assert loc is None or frozenset(loc) in cand
