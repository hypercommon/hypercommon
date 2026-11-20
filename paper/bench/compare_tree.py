"""
The per-branch threshold (hc_tree, results_tree/) against HyperCommon and the other
methods (results2/), exploratory. Protocols and noise-averaged Omega as analysis2.py.

    python compare_tree.py [experiments...]
Writes results_tree/compare_bench.csv and prints one row per family and protocol.
"""
import os
import sys

import numpy as np
import pandas as pd

import analysis2 as A

HERE = os.path.dirname(os.path.abspath(__file__))
EXPS = ["rings_uniform", "rings_hetero", "rings_hetero_overlap", "rings_pref",
        "lfr", "lfr_slices", "lfr_clustered"]
ABL = {"hc_2test", "hc_tri", "hc_econn"}


def load(exp):
    parts = []
    for f in (f"results2/{exp}.csv", f"results2/{exp}_oslom.csv", f"results2/{exp}_sbm.csv",
              f"results_tree/{exp}.csv"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            parts.append(pd.read_csv(p, dtype={"param": str}, on_bad_lines="skip"))
    d = pd.concat(parts, ignore_index=True)
    d = d[d.error.isna() & (d.split == "test")].copy()
    d["param"] = d["param"].fillna("")
    return d


def main(exps):
    rows = []
    for exp in exps:
        d = load(exp)
        d = d[~d.method.isin(A.NATIVE_OF)]
        key = ["family", "method", "p", "seed"]
        orc = d.loc[d.groupby(key).omega.idxmax()].assign(protocol="oracle")
        lf = d.loc[d.groupby(key).eq.idxmax()].assign(protocol="label_free")
        # the tree at one fixed minimum size, no selection
        fix = d[(d.method == "hc_tree") & (d.param == "50")].assign(protocol="fixed50")
        sel = pd.concat([orc, lf], ignore_index=True)
        for (fam, prot), g in sel.groupby(["family", "protocol"]):
            area = {m: A.areas(gm, "omega") for m, gm in g.groupby("method")}
            area = {m: a for m, a in area.items() if a is not None}
            if "hc_tree" not in area:
                continue
            mean = pd.Series({m: a.mean() for m, a in area.items()})
            others = mean.drop([m for m in mean.index if m in ABL or m in ("hc", "hc_tree")])
            ranked = mean.drop([m for m in mean.index if m in ABL or m == "hc"]).sort_values(ascending=False)
            bo = others.idxmax()
            dd = (area["hc_tree"] - area[bo]).dropna()
            dh = (area["hc_tree"] - area["hc"]).dropna()
            f50 = A.areas(fix[fix.family == fam], "omega")
            rows.append(dict(experiment=exp, family=fam, protocol=prot, tree=mean["hc_tree"],
                             hc=mean["hc"], best_other=bo, best_other_mean=others[bo],
                             tree_minus_hc=dh.mean(), tree_minus_best=dd.mean(),
                             tree_rank=list(ranked.index).index("hc_tree") + 1,
                             hc_rank=list(mean.drop([m for m in mean.index if m in ABL or m == "hc_tree"])
                                          .sort_values(ascending=False).index).index("hc") + 1,
                             tree_fixed50=f50.mean() if f50 is not None else np.nan))
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(HERE, "results_tree", "compare_bench.csv"), index=False)
    with pd.option_context("display.width", 250):
        print(r.round(3).to_string(index=False))
        for prot, g in r.groupby("protocol"):
            ring = g[~g.experiment.str.startswith("lfr")]
            lfr = g[g.experiment.str.startswith("lfr")]
            for name, x in (("rings", ring), ("lfr", lfr)):
                print(f"{prot:10s} {name:5s} n={len(x)}  tree>hc {int((x.tree_minus_hc > 0.005).sum())}  "
                      f"tree<hc {int((x.tree_minus_hc < -0.005).sum())}  mean tree-hc {x.tree_minus_hc.mean():+.3f}  "
                      f"tree rank1 {int((x.tree_rank == 1).sum())}  hc rank1 {int((x.hc_rank == 1).sum())}  "
                      f"tree top3 {int((x.tree_rank <= 3).sum())}  hc top3 {int((x.hc_rank <= 3).sum())}")


if __name__ == "__main__":
    main(sys.argv[1:] or EXPS)
