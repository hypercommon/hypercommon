"""Robustness to noise on the ring families: per family and method, under the oracle
and label-free protocols, the mean Omega over test graphs at each p, its mean over the
high-noise half (p >= 0.4), and p_half, the first p at which the mean drops below 0.5
(linear interpolation; '>0.8' if never). Writes results2/robust_p.csv."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results2")
EXPS = ["rings_uniform", "rings_hetero", "rings_hetero_overlap", "rings_pref"]


def p_half(ps, ys, level=0.5):
    for i in range(len(ps)):
        if ys[i] < level:
            if i == 0:
                return ps[0]
            a, b = ys[i - 1], ys[i]
            return ps[i - 1] + (a - level) / (a - b) * (ps[i] - ps[i - 1])
    return np.inf


out = []
for e in EXPS:
    d = pd.read_csv(os.path.join(R, f"{e}.csv"), usecols=["family", "p", "seed", "split", "method", "omega", "eq"])
    d = d[(d.split == "test") & d.omega.notna()]
    key = ["family", "method", "p", "seed"]
    for prot, col in (("oracle", "omega"), ("label_free", "eq")):
        sel = d.loc[d.groupby(key)[col].idxmax()]
        m = sel.groupby(["family", "method", "p"]).omega.mean().reset_index()
        for (fam, meth), g in m.groupby(["family", "method"]):
            g = g.sort_values("p")
            ps, ys = g.p.values, g.omega.values
            hi = g[g.p >= 0.4]
            out.append(dict(experiment=e, family=fam, method=meth, protocol=prot,
                            high_noise_mean=float(np.trapezoid(hi.omega, hi.p) / (hi.p.max() - hi.p.min())),
                            omega_at_max_p=float(ys[-1]), p_half=p_half(ps, ys)))
pd.DataFrame(out).to_csv(os.path.join(R, "robust_p.csv"), index=False)
