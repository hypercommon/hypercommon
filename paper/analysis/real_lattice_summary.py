"""
Summary of the detection runs on the lattice-like real networks
(bench/real_lattice_detection.py): Hi-C contact maps with contact domains and protein
contact networks with CATH domains, each with its distance- or sequence-only control.

For every graph and method: the best value over the parameter grid (best) and the value at
the parameter of largest overlapping modularity (free), stochastic methods averaged over
seeds, for Omega and F1. Per set of graphs: the mean of each method, and HyperCommon minus
the best other method chosen per set (the method with the highest mean), with a 95%
bootstrap interval over the graphs (paired by graph).

    python real_lattice_summary.py [kinds...]     (default: hic protein)
Writes figures/data/real_lattice_summary.csv and prints the tables.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "bench", "results")
OUT = os.path.join(HERE, "..", "figures", "data", "real_lattice_summary.csv")
NOT_OTHER = {"hc", "hc_tri", "hc_2test", "hc_econn", "components"}
rng = np.random.default_rng(0)


def load(kind):
    parts = [pd.read_csv(os.path.join(RES, f), dtype={"param": str}, low_memory=False)
             for f in (f"real_lattice_{kind}.csv", f"real_lattice_{kind}_wsl.csv")
             if os.path.exists(os.path.join(RES, f)) and os.path.getsize(os.path.join(RES, f)) > 0]
    d = pd.concat(parts)
    return d[d.error.isna()]


def picks(d, metric):
    key = ["network", "method", "seed"]
    best = d.loc[d.groupby(key)[metric].idxmax()].groupby(["network", "method"])[metric].mean()
    free = d.loc[d.groupby(key)["eq"].idxmax()].groupby(["network", "method"])[metric].mean()
    return best.unstack(), free.unstack()


def setname(kind, net):
    if kind == "hic":
        k = net.split("_k")[1].split("_")[0]
        return f"Hi-C k={k}" + (" control" if net.endswith("_null") else "")
    if kind == "protein":
        return "protein " + ("7A" if net.endswith("_7A") else "control" if net.endswith("_null") else "8A")
    base = net[:-5] if net.endswith("_null") else net
    parts = base.split("_")
    k = [p[1:] for p in parts if p[:1] == "k" and p[1:].isdigit()]
    if kind == "geo":                            # each geo graph is its own set
        tag = [p for p in parts if not (p[:1] == "k" and p[1:].isdigit())]
    else:                                        # graph type words (tssb knn/band, spatial parts)
        tag = [p for p in parts if p in ("knn", "band", "parts") or (p[:1] == "r" and p[1:].isdigit())]
    return " ".join([kind] + tag) + (f" k={k[-1]}" if k else "") + (" control" if net.endswith("_null") else "")


def main():
    rows = []
    import sys
    for kind in (sys.argv[1:] or ("hic", "protein")):
        try:
            d = load(kind)
        except ValueError:
            continue
        d = d[~d.method.str.startswith("hcf_")]
        for metric in ("omega", "f1"):
            for prot, T in zip(("best", "free"), picks(d, metric)):
                T = T.copy()
                T["set"] = [setname(kind, n) for n in T.index]
                for s, g in T.groupby("set"):
                    g = g.drop(columns="set").dropna(axis=1, how="all")
                    others = [m for m in g.columns if m not in NOT_OTHER]
                    means = g.mean()
                    bo = means[others].idxmax()
                    diff = (g["hc"] - g[bo]).dropna().to_numpy()
                    bs = [rng.choice(diff, len(diff)).mean() for _ in range(2000)]
                    lo, hi = np.percentile(bs, [2.5, 97.5])
                    wins = int((g["hc"] > g[others].max(axis=1) + 1e-9).sum())
                    rows.append(dict(set=s, metric=metric, protocol=prot, graphs=len(g), hc=means["hc"],
                                     tri=means.get("hc_tri"), two_test=means.get("hc_2test"),
                                     best_other=bo, best_other_mean=means[bo], diff=diff.mean(), lo=lo, hi=hi,
                                     hc_best_on=wins,
                                     ranking=", ".join(f"{m} {v:.2f}" for m, v in
                                                       means.drop(["components"], errors="ignore")
                                                       .sort_values(ascending=False).head(6).items())))
    out = pd.DataFrame(rows)
    out.to_csv(OUT if not sys.argv[1:] else OUT.replace(".csv", "_" + "_".join(sys.argv[1:]) + ".csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 90):
        print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
