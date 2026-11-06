import random, sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import networkx as nx
from hc_fast import TripleScores
from hypercommon.algorithm import get_communities
from predicates.jaccard import closed_neighborhood_jaccard_predicate as P

def canon(cs): return sorted(sorted(c) for c in cs)

def test_equivalence():
    rng = random.Random(0)
    graphs = [nx.karate_club_graph(), nx.les_miserables_graph()]
    graphs = [nx.convert_node_labels_to_integers(g) for g in graphs]
    for s in range(6):
        graphs.append(nx.watts_strogatz_graph(150, rng.choice([4, 6, 8]), 0.2, seed=s))
        graphs.append(nx.gnp_random_graph(80, 0.08, seed=s))
    for G in graphs:
        ts = TripleScores(G)
        for t in [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0]:
            assert canon(ts.communities("hc", t)) == canon(get_communities(G, P(t))), t

def test_triangles_only_is_filtered_cpm():
    from networkx.algorithms.community import k_clique_communities
    from predicates.jaccard import closed_neighborhood_jaccard
    from hypercommon.hypernode import HCNode
    for s in range(5):
        G = nx.watts_strogatz_graph(120, 8, 0.3, seed=s)
        hc = {u: HCNode(u, set(G[u])) for u in G}
        ts = TripleScores(G)
        for t in [0.0, 0.1, 0.2, 0.3]:
            F = nx.Graph([(u, v) for u, v in G.edges() if closed_neighborhood_jaccard(hc[u], hc[v]) >= t])
            ref = [set(c) for c in k_clique_communities(F, 3)] if F.number_of_edges() else []
            assert canon(ts.communities("hc_tri", t)) == canon(ref), (s, t)
