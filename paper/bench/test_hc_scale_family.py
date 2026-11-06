import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import networkx as nx
import numpy as np

from hc_family import FUNCTIONS, FamilyScores
from hc_scale_family import ScaleFamily


def canon(cs):
    return sorted(sorted(c) for c in cs)


def test_family_scale_equals_reference():
    graphs = [nx.karate_club_graph()]
    for s in range(4):
        graphs.append(nx.watts_strogatz_graph(150, 8, 0.25, seed=s))
        graphs.append(nx.barabasi_albert_graph(120, 3, seed=s))
    for G in graphs:
        a, b = FamilyScores(G), ScaleFamily(G)
        for f in FUNCTIONS:
            # same multiset of triple scores
            sa = np.sort(a.triple_score(f))
            sb = np.sort(b.triple_score(f))
            assert np.allclose(sa, sb), f
            for q in (0.0, 0.3, 0.6, 0.9):
                assert canon(a.communities_f(f, q)) == canon(b.communities_f(f, q)), (f, q)
