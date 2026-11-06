import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import networkx as nx

from hc_fast import TripleScores
from hc_scale import ScaleScores


def canon(cs):
    return sorted(sorted(c) for c in cs)


def graphs():
    out = [nx.karate_club_graph(), nx.convert_node_labels_to_integers(nx.les_miserables_graph())]
    for s in range(5):
        out.append(nx.watts_strogatz_graph(200, 8, 0.25, seed=s))
        out.append(nx.gnp_random_graph(120, 0.06, seed=s))
        out.append(nx.barabasi_albert_graph(150, 3, seed=s))
    return out


def test_same_triples_and_scores():
    for G in graphs():
        a, b = TripleScores(G), ScaleScores(G)
        ka = {tuple(sorted(G_ids)): s for G_ids, s in zip(
            [tuple(a.nodes[x] for x in t) for t in a.triples], a.s3)}
        kb = {tuple(sorted(b.nodes[x] for x in t)): float(s) for t, s in zip(b.triples, b.s3)}
        assert ka.keys() == kb.keys()
        assert max(abs(ka[k] - kb[k]) for k in ka) < 1e-6
        assert a.n_wedges == b.n_wedges


def test_same_communities():
    for G in graphs():
        a, b = TripleScores(G), ScaleScores(G)
        for t in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]:
            assert canon(b.communities(t)) == canon(a.communities("hc", t)), t
            assert canon(b.communities(t, triangles_only=True)) == canon(a.communities("hc_tri", t)), t
