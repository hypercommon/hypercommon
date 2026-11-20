"""
The per-branch threshold (hc_tree, results_tree/real_<kind>.csv) on the real networks
against HyperCommon and the fourteen other methods, exploratory. Picks as
analysis/real_lattice_summary.py: oracle = best parameter per graph, label-free = largest
overlapping modularity, stochastic methods averaged over seeds.

    python compare_tree_real.py
Writes results_tree/compare_real.csv.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import real_lattice_summary as S  # noqa: E402

KINDS = ["hic", "sketch", "image", "tracks", "geo", "protein", "transit", "streets", "tssb",
         "spatial", "text"]
NOT_RANKED = {"hc_tri", "hc_2test", "hc_econn", "components"}


def tree(kind):
    p = os.path.join(HERE, "results_tree", f"real_{kind}.csv")
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p, dtype={"param": str})
    return d[d.error.isna()]


def compare(d, setname):
    rows = []
    d = d[~d.method.str.startswith("hcf_") & ~d.method.isin(NOT_RANKED)]
    for metric in ("omega", "f1"):
        for prot, T in zip(("oracle", "label_free"), S.picks(d, metric)):
            T = T.copy()
            T["set"] = [setname(n) for n in T.index]
            for s, g in T.groupby("set"):
                g = g.drop(columns="set").dropna(axis=1, how="all")
                if "hc" not in g or "hc_tree" not in g:
                    continue
                g = g.dropna(subset=["hc", "hc_tree"])
                if g.empty:
                    continue
                m = g.mean()
                others = m.drop(["hc", "hc_tree"])
                bo = others.idxmax()
                rank_t = list(m.drop("hc").sort_values(ascending=False).index).index("hc_tree") + 1
                rank_h = list(m.drop("hc_tree").sort_values(ascending=False).index).index("hc") + 1
                rows.append(dict(set=s, metric=metric, protocol=prot, graphs=len(g),
                                 tree=m["hc_tree"], hc=m["hc"], best_other=bo, best_other_mean=others[bo],
                                 tree_minus_hc=m["hc_tree"] - m["hc"], tree_rank=rank_t, hc_rank=rank_h,
                                 methods=len(m) - 1))
    return rows


def main():
    rows = []
    for kind in KINDS:
        t = tree(kind)
        if t is None:
            continue
        d = pd.concat([S.load(kind), t], ignore_index=True)
        d = d[~d.network.str.endswith("_null")]        # real graphs only
        for r in compare(d, lambda n, k=kind: S.setname(k, n)):
            rows.append(dict(kind=kind, **r))
    t = tree("more")
    if t is not None:
        base = pd.read_csv(os.path.join(HERE, "results", "real_more_v3.csv"), dtype={"param": str},
                           low_memory=False)
        d = pd.concat([base[base.error.isna()], t], ignore_index=True)
        for r in compare(d, lambda n: n):
            rows.append(dict(kind="more", **r))
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(HERE, "results_tree", "compare_real.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(r[r.metric == "omega"].round(3).to_string(index=False))
        for (metric, prot), g in r.groupby(["metric", "protocol"]):
            for name, x in (("lattice", g[g.kind != "more"]), ("more", g[g.kind == "more"])):
                print(f"{metric:5s} {prot:10s} {name:8s} sets={len(x)}  tree>hc {int((x.tree_minus_hc > 0.005).sum())} "
                      f"tree<hc {int((x.tree_minus_hc < -0.005).sum())}  mean {x.tree_minus_hc.mean():+.3f}  "
                      f"median rank tree {x.tree_rank.median():g} hc {x.hc_rank.median():g}  "
                      f"rank1 tree {int((x.tree_rank == 1).sum())} hc {int((x.hc_rank == 1).sum())}  "
                      f"top3 tree {int((x.tree_rank <= 3).sum())} hc {int((x.hc_rank <= 3).sum())}")


if __name__ == "__main__":
    main()
