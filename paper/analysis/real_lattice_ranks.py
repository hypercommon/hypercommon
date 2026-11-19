"""
HyperCommon's rank among the fifteen methods (HyperCommon and the fourteen others; the
ablations are left out) on the lattice-like real networks of real_lattice_detection.py.

For every set of graphs (a data kind at one k, or one geographic graph at one k), every
method's mean Omega and F1 over the set's graphs, under the oracle (best parameter per
graph) and the label-free (largest overlapping modularity) protocols, as in
real_lattice_summary.py; HyperCommon's rank among the fifteen means (1 = most accurate),
the most accurate other method and HyperCommon minus it with a paired 95% bootstrap
interval over the graphs. Controls (the graphs with the data signal removed) are ranked
too and marked.

    python real_lattice_ranks.py
Writes figures/data/real_lattice_ranks.csv, sections/table_real_lattice.tex (one row per
kind, main text) and sections/table_real_lattice_sets.tex (one row per set, supplement).
"""
import os

import numpy as np
import pandas as pd

from real_lattice_summary import load, picks, setname

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "data", "real_lattice_ranks.csv")
TEX = os.path.join(HERE, "..", "sections", "table_real_lattice.tex")
TEX_SETS = os.path.join(HERE, "..", "sections", "table_real_lattice_sets.tex")
NOT_RANKED = {"hc_tri", "hc_2test", "hc_econn", "components"}

# kind, name in the tables, what a group is; in the order of the tables
KINDS = [
    ("hic", "Hi-C maps", "domain"),
    ("sketch", "sketches", "stroke"),
    ("image", "images", "segment"),
    ("tracks", "pedestrians", "track"),
    ("geo", "commuting", "region"),
    ("protein", "proteins", "domain"),
    ("transit", "transit", "route"),
    ("streets", "streets", "street"),
    ("tssb", "time series", "segment"),
    ("spatial", "brain tissue", "layer"),
    ("text", "articles", "section"),
]
SHORT = {"leiden": "Leiden", "leiden_cpm": "L-CPM", "linkcomm": "Link", "cpm": "CPM",
         "scan": "SCAN", "scan_j": "SCAN-J", "infomap": "Infomap", "walktrap": "Walktrap",
         "lpa": "LPA", "slpa": "SLPA", "demon": "DEMON", "angel": "ANGEL", "oslom": "OSLOM",
         "sbm": "SBM", "hc": "HC"}
rng = np.random.default_rng(0)


def ranks():
    rows = []
    for kind, _, _ in KINDS:
        d = load(kind)
        d = d[~d.method.str.startswith("hcf_")]
        base = d.network[~d.network.str.endswith("_null")].str.replace(r"(_(k\d+|7A|r\d|knn|band|parts))+$", "", regex=True)
        n_base = base.nunique()
        for metric in ("omega", "f1"):
            for prot, T in zip(("oracle", "label_free"), picks(d, metric)):
                T = T[[c for c in T.columns if c not in NOT_RANKED]].copy()
                T["set"] = [setname(kind, n) for n in T.index]
                for s, g in T.groupby("set"):
                    g = g.drop(columns="set").dropna(axis=1, how="all")
                    m = g.mean().sort_values(ascending=False)
                    others = m.drop("hc")
                    bo = others.index[0]
                    diff = (g["hc"] - g[bo]).dropna().to_numpy()
                    bs = rng.choice(diff, (2000, len(diff))).mean(1)
                    lo, hi = np.percentile(bs, [2.5, 97.5])
                    rows.append(dict(kind=kind, set=s, control=s.endswith("control"), metric=metric,
                                     protocol=prot, graphs=len(g), networks=n_base, methods=len(m),
                                     rank=list(m.index).index("hc") + 1, hc=m["hc"],
                                     best_other=bo, best_other_mean=others.iloc[0],
                                     diff=diff.mean(), lo=lo, hi=hi))
    return pd.DataFrame(rows)


def _k(s):
    """Sort key: the set's tag, then k."""
    base = s.replace(" control", "")
    k = int(base.split("k=")[1]) if "k=" in base else 0
    return base.split("k=")[0].strip(), k


def table_kinds(r):
    real = r[~r.control]
    lines = [
        r"\begin{table}[t]",
        r"\caption{HyperCommon (HC) on lattice-like real networks: its rank by $\Omega$ among "
        r"the fifteen methods (HC and the fourteen others; 1 is the most accurate), with the "
        r"best parameter per graph (oracle) and without labels (label-free). A kind holds one "
        r"or more sets of graphs (one per $k$, and per graph type or map); a single number is "
        r"the rank on every set, two numbers the best and the median rank over the sets. "
        r"Graphs: distinct networks of the kind. Leader: the most accurate method with the oracle on most sets (several when tied). Every set, F1 and "
        r"the controls are in Table~\ref{tab:S-real-lattice-sets}.}",
        r"\label{tab:real-lattice}",
        r"\scriptsize\setlength{\tabcolsep}{3pt}",
        r"\begin{ruledtabular}",
        r"\begin{tabular}{llrccl}",
        r"Kind & group & graphs & oracle & free & leader \\",
        r"\hline",
    ]
    for kind, name, group in KINDS:
        o = real[(real.kind == kind) & (real.metric == "omega")]
        cells = []
        for prot in ("oracle", "label_free"):
            x = o[o.protocol == prot]["rank"]
            cells.append(f"{x.min()}" if x.min() == x.max() else f"{x.min()}, {x.median():g}")
        orc = o[o.protocol == "oracle"]
        vc = orc.apply(lambda q: "hc" if q["rank"] == 1 else q.best_other, axis=1).value_counts()
        lead = "/".join(SHORT[m] for m in vc.index[vc == vc.iloc[0]])  # ties: all
        lines.append(f"{name} & {group} & {int(orc.networks.iloc[0])} & {cells[0]} & "
                     f"{cells[1]} & {lead} \\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    return "\n".join(lines) + "\n"


def table_sets(r):
    lines = [
        r"\begin{longtable*}{llrrrrrlr}",
        r"\caption{Every set of lattice-like real graphs: HyperCommon's rank among the fifteen "
        r"methods by $\Omega$ and F1, oracle (orc) and label-free (free) protocols; its mean "
        r"$\Omega$ (oracle), the most accurate other method with the oracle and HC minus it, "
        r"with the paired 95\% bootstrap interval over the graphs of the set. Controls keep "
        r"the nodes and groups and build the edges from position alone (sequence, genomic, "
        r"grid or map distance; time for sketches), removing the data signal.}"
        r"\label{tab:S-real-lattice-sets}\\",
        r"\hline\hline",
        r"Kind & set & $\Omega$ orc & $\Omega$ free & F1 orc & F1 free & HC $\Omega$ & best other & "
        r"difference \\",
        r"\hline\endfirsthead",
        r"\hline Kind & set & $\Omega$ orc & $\Omega$ free & F1 orc & F1 free & HC $\Omega$ & best other "
        r"& difference \\ \hline\endhead",
        r"\hline\hline\endlastfoot",
    ]
    for kind, name, _ in KINDS:
        sub = r[r.kind == kind]

        def tag(s):
            t = s.replace(" control", "").replace("Hi-C", "").replace(kind, "", 1).strip()
            return t or "8A"                       # the protein control is built on the 8 A graphs

        first = True
        for s in sorted(sub.set.unique(), key=lambda s: (_k(tag(s)), s.endswith("control"))):
            q = sub[sub.set == s].set_index(["metric", "protocol"])
            ro = q.loc[("omega", "oracle")]
            ctl = "control" if s.endswith("control") else ""
            lines.append(
                f"{name if first else ''} & {tag(s) + (', control' if ctl else '')} & "
                f"{q.loc[('omega', 'oracle'), 'rank']} & {q.loc[('omega', 'label_free'), 'rank']} & "
                f"{q.loc[('f1', 'oracle'), 'rank']} & {q.loc[('f1', 'label_free'), 'rank']} & "
                f"{ro.hc:.2f} & {SHORT[ro.best_other]} {ro.best_other_mean:.2f} & "
                f"${ro['diff']:+.2f}$ [${ro.lo:+.2f}, {ro.hi:+.2f}$] \\\\")
            first = False
        lines.append(r"\hline")
    lines.append(r"\end{longtable*}")
    return "\n".join(lines) + "\n"


def main():
    r = ranks()
    r.to_csv(OUT, index=False)
    with open(TEX, "w", newline="\n") as fh:
        fh.write(table_kinds(r))
    with open(TEX_SETS, "w", newline="\n") as fh:
        fh.write(table_sets(r))
    real = r[~r.control & (r.metric == "omega")]
    w = real.pivot_table(index="set", columns="protocol", values="rank")
    print(w.to_string())
    top3 = real[(real.protocol == "oracle")].groupby("kind")["rank"].min()
    print("kinds with HC in the top 3 on some set (oracle, Omega):", (top3 <= 3).sum(), "of", len(top3))


if __name__ == "__main__":
    main()
