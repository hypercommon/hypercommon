"""
Table of the per-branch threshold against HyperCommon (Supplemental Material; the "first
try" of the Discussion), from results_tree/compare_bench.csv (compare_tree.py) and
results_tree/compare_real.csv (compare_tree_real.py).

    python make_table_tree.py
Writes ../sections/table_tree.tex.
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "results_tree")


def main():
    b = pd.read_csv(os.path.join(D, "compare_bench.csv"))
    b["grp"] = b.experiment.str.startswith("lfr").map({True: "lfr", False: "rings"})
    r = pd.read_csv(os.path.join(D, "compare_real.csv"))
    r = r[r.metric == "omega"].copy()
    r["grp"] = (r.kind == "more").map({True: "more", False: "lattice"})
    a = pd.concat([b[["grp", "protocol", "tree", "hc", "tree_minus_hc", "tree_rank", "hc_rank"]],
                   r[["grp", "protocol", "tree", "hc", "tree_minus_hc", "tree_rank", "hc_rank"]]])
    rows = [("rings", "Ring lattices (settings)"), ("lfr", "LFR (families)"),
            ("lattice", "Lattice-like real (sets)"), ("more", "Fourteen networks")]
    lines = []
    for g, name in rows:
        cells = []
        for prot in ("oracle", "label_free"):
            x = a[(a.grp == g) & (a.protocol == prot)]
            up, dn = int((x.tree_minus_hc > 0.005).sum()), int((x.tree_minus_hc < -0.005).sum())
            cells.append(f"{x.hc.mean():.3f} & {x.tree.mean():.3f} & {up}/{dn} & "
                         f"{x.hc_rank.median():g} & {x.tree_rank.median():g}")
        n = len(a[(a.grp == g) & (a.protocol == "oracle")])
        lines.append(f"{name}, {n} & " + " & ".join(cells) + " \\\\")

    def cell(fam, prot):
        x = b[(b.family == fam) & (b.protocol == prot)].iloc[0]
        return x.hc, x.tree, x.best_other_mean
    g0 = cell("gradient10_o00", "label_free")
    g5 = cell("gradient10", "label_free")
    cap = ("The per-branch threshold (\\texttt{bench/hc\\_tree.py}) against HyperCommon (HC), $\\Omega$. "
           "The covers of HyperCommon at all thresholds form one tree (Kruskal on the triples in "
           "decreasing $s_3$, over the pairs); it is condensed with a minimum cluster size, in "
           "pairs, and clusters are selected by excess of mass with their stability counted in "
           "nodes, as in HDBSCAN~\\cite{campello2013density}. The minimum size is the parameter "
           "(7 values, 5--500), chosen as for every method: best per graph (oracle) or largest "
           "overlapping modularity (label-free). Rings and LFR: mean over settings of the "
           "noise-averaged $\\Omega$; real networks: mean over sets (lattice-like kinds, "
           "Table~\\ref{tab:S-real-lattice-sets}) or networks. Up/down: settings, sets or networks "
           "where the per-branch threshold is above / below HC by more than 0.005. Rank: median "
           "rank among HC or the per-branch threshold and the other methods. Its oracle chooses "
           "among 7 values, HC's among 30 on the benchmarks and 100 on the real networks. On the "
           f"gradient without labels it lifts HC from {g0[0]:.2f} to {g0[1]:.2f} without overlap "
           f"and from {g5[0]:.2f} to {g5[1]:.2f} at $o_n=0.05$ (link communities: {g0[2]:.2f} "
           f"and {g5[2]:.2f}).")
    tex = "\n".join([
        "\\begin{table*}[t]", "\\caption{" + cap + "}", "\\label{tab:S-tree}",
        "\\begin{ruledtabular}", "\\begin{tabular}{lccccc|ccccc}",
        " & \\multicolumn{5}{c}{Oracle} & \\multicolumn{5}{c}{Label-free} \\\\",
        " & HC & Branch & Up/down & Rank HC & Rank branch & HC & Branch & Up/down & Rank HC & Rank branch \\\\",
        "\\hline", *lines, "\\end{tabular}", "\\end{ruledtabular}", "\\end{table*}"])
    with open(os.path.join(HERE, "..", "sections", "table_tree.tex"), "w") as fh:
        fh.write(tex + "\n")
    print(tex)


if __name__ == "__main__":
    main()
