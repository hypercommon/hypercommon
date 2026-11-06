"""
Why HyperCommon loses on LFR: is the pairwise signal weak, or does the construction
fail to use it? For LFR test graphs (standard o_n = 0.1, and the slice with larger
communities, sizes 50-250), at every mixing level, seeds 0-4:

  auc         AUC of closed-neighbourhood Jaccard for edges inside a planted community
              against edges between communities (the pairwise signal)
  pure_adm_*  share of the wedges lying inside one planted community that pass the
              three-pair test at t = 0.1, for communities of size <= 40 and >= 80
              (whether large, internally sparse communities can stay connected)
  mix_adm     share of mixed wedges (spanning two communities) passing at t = 0.1

Writes figures/data/lfr_signal.csv and prints family x mu means.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench"))
import experiments  # noqa: E402
from hc_fast import TripleScores  # noqa: E402

T = 0.1


def auc(pos, neg):
    from scipy.stats import mannwhitneyu
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


def main():
    rows = []
    for exp, fam in (("lfr", "on10"), ("lfr_slices", "big")):
        E = experiments.EXPERIMENTS[exp]
        for mu in experiments.MU_GRID:
            for seed in range(5):
                G, truth = E["build"](E["families"][fam], mu, seed)
                memb = {}
                for i, c in enumerate(truth):
                    for u in c:
                        memb.setdefault(u, set()).add(i)
                N = {u: set(G[u]) | {u} for u in G}
                inner, outer = [], []
                for u, v in G.edges():
                    j = len(N[u] & N[v]) / len(N[u] | N[v])
                    (inner if memb.get(u, set()) & memb.get(v, set()) else outer).append(j)
                ts = TripleScores(G)
                nodes = ts.nodes
                com_of = [memb.get(u, set()) for u in nodes]
                size = {i: len(c) for i, c in enumerate(truth)}
                shared = [com_of[a] & com_of[b] & com_of[c] for a, b, c in ts.triples]
                adm = ts.s3 >= T
                small = np.array([bool(s) and min(size[i] for i in s) <= 40 for s in shared])
                large = np.array([bool(s) and min(size[i] for i in s) >= 80 for s in shared])
                mixed = np.array([not s for s in shared])
                rows.append(dict(family=f"{exp}:{fam}", mu=mu, seed=seed, auc=auc(inner, outer),
                                 pure_adm_small=adm[small].mean() if small.any() else np.nan,
                                 pure_adm_large=adm[large].mean() if large.any() else np.nan,
                                 mix_adm=adm[mixed].mean() if mixed.any() else np.nan))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(HERE, "..", "figures", "data", "lfr_signal.csv"), index=False)
    print(d.groupby(["family", "mu"]).mean(numeric_only=True).drop(columns="seed").round(3).to_string())


if __name__ == "__main__":
    main()
