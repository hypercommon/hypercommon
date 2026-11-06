"""
Where the triangles of LFR graphs lie: the fraction of triangles whose three nodes share
a planted community, for the standard overlapping LFR family (o_n = 0.1) and the two
clustered ones (generator option C = 0.3, 0.5), at every mixing level, test seeds 0-9.
Writes figures/data/lfr_triangles.csv (one row per graph) and prints family means.
"""
import os
import sys

import networkx as nx
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(HERE, "..", "bench"))
import experiments  # noqa: E402

CASES = [("lfr", "on10"), ("lfr_clustered", "on10_C03"), ("lfr_clustered", "on10_C05")]


def triangle_stats(G, truth):
    memb = {}
    for i, c in enumerate(truth):
        for u in c:
            memb.setdefault(u, set()).add(i)
    inside = total = 0
    for u in G:
        nu = set(G[u])
        for v in nu:
            if v <= u:
                continue
            for w in nu & set(G[v]):
                if w <= v:
                    continue
                total += 1
                inside += bool(memb.get(u, set()) & memb.get(v, set()) & memb.get(w, set()))
    return total, inside


def main():
    rows = []
    for exp, fam in CASES:
        E = experiments.EXPERIMENTS[exp]
        for mu in E["p_grid"]:
            for seed in range(10):
                G, truth = E["build"](E["families"][fam], mu, seed)
                tot, ins = triangle_stats(G, truth)
                rows.append(dict(experiment=exp, family=fam, mu=mu, seed=seed, triangles=tot,
                                 inside=ins, frac_inside=ins / tot if tot else float("nan"),
                                 transitivity=nx.transitivity(G)))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(HERE, "..", "figures", "data", "lfr_triangles.csv"), index=False)
    print(d.groupby(["family", "mu"])[["triangles", "frac_inside", "transitivity"]].mean().round(3).to_string())


if __name__ == "__main__":
    main()
