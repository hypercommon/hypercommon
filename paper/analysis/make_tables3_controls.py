"""
Control families: rings with hubs, split overlap, the
degree-preserving geometry control, random geometric communities, and the null graphs.

  sections/supp_controls.tex     tab:S-controls (Omega, oracle and label-free, with the
                                 reference family each control changes) and tab:S-nulls
  figures/data/controls3.csv     the numbers, for the text
Sources: bench/results2/<exp>_summary2.csv and _paired2.csv (analysis2.py).

    python make_tables3_controls.py
"""
import os

import numpy as np
import pandas as pd

from make_tables3 import NAME, OTHER_IDS, R2, SEC, DATA

# (experiment, family, label, reference experiment, reference family, reference label)
CONTROLS = [
    ("rings_hubs", "hub_s50", r"hubs, $(50,k_i\in[4,32])\times40$", "rings_uniform", "k8_s50_o05", r"uniform $k{=}8$"),
    ("rings_hubs", "hub_s100", r"hubs, $(100,k_i\in[8,64])\times20$", "rings_uniform", "k16_s100_o05", r"uniform $k{=}16$"),
    ("rings_hubs", "hub_two_worlds", "hubs, two worlds", "rings_hetero", "two_worlds", "two worlds"),
    ("rings_split", "two_worlds_o05_split", r"split overlap, two worlds, $o_n{=}0.05$", "rings_hetero", "two_worlds", "two worlds"),
    ("rings_split", "two_worlds_o10_split", r"split overlap, two worlds, $o_n{=}0.10$", "rings_hetero_overlap", "two_worlds_o10", r"two worlds, $o_n{=}0.10$"),
    ("rings_split", "three_big_many_small_o05_split", r"split overlap, three big, $o_n{=}0.05$", "rings_hetero", "three_big_many_small", "three big"),
    ("rings_split", "three_big_many_small_o10_split", r"split overlap, three big, $o_n{=}0.10$", "rings_hetero_overlap", "three_big_many_small_o10", r"three big, $o_n{=}0.10$"),
    ("rings_geometry_swap", "k8_s50_o05_q010", r"swaps, uniform $k{=}8$, $q{=}0.1$", "rings_geometry", "k8_s50_o05_q010", "endpoint, $q{=}0.1$"),
    ("rings_geometry_swap", "k8_s50_o05_q030", r"swaps, uniform $k{=}8$, $q{=}0.3$", "rings_geometry", "k8_s50_o05_q030", "endpoint, $q{=}0.3$"),
    ("rings_geometry_swap", "three_big_many_small_q010", r"swaps, three big, $q{=}0.1$", "rings_geometry", "three_big_many_small_q010", "endpoint, $q{=}0.1$"),
    ("rings_geometry_swap", "three_big_many_small_q030", r"swaps, three big, $q{=}0.3$", "rings_geometry", "three_big_many_small_q030", "endpoint, $q{=}0.3$"),
    ("rgg", "rgg_uniform", r"random geometric, $(50,8)\times40$", "rings_uniform", "k8_s50_o05", r"uniform $k{=}8$"),
    ("rgg", "rgg_two_worlds", "random geometric, two worlds", "rings_hetero", "two_worlds", "two worlds"),
]
CTRL_OTHERS = [m for m in OTHER_IDS if m not in ("demon", "angel", "scan", "scan_j")]


def load(exps):
    S, P = [], []
    for e in sorted(set(exps)):
        f = os.path.join(R2, f"{e}_summary2.csv")
        if os.path.exists(f):
            S.append(pd.read_csv(f).assign(experiment=e))
            P.append(pd.read_csv(os.path.join(R2, f"{e}_paired2.csv")).assign(experiment=e))
    return pd.concat(S, ignore_index=True), pd.concat(P, ignore_index=True)


def cell(S, P, exp, fam, prot, p_max=None):
    g = S[(S.experiment == exp) & (S.family == fam) & (S.protocol == prot) & (S.metric == "omega")]
    g = g.set_index("method")["mean"]
    if g.empty or "hc" not in g:
        return None
    o = g.reindex(CTRL_OTHERS).dropna()
    bm = o.idxmax()
    pr = P[(P.experiment == exp) & (P.family == fam) & (P.protocol == prot) & (P.metric == "omega") & (P.method == bm)]
    mark = "=" if pr.empty else ("+" if pr.hi.iloc[0] < 0 else "-" if pr.lo.iloc[0] > 0 else "=")
    return dict(hc=g["hc"], tri=g.get("hc_tri", np.nan), best=o[bm], best_m=bm, mark=mark)


def main():
    exps = [c[0] for c in CONTROLS] + [c[3] for c in CONTROLS] + ["nulls"]
    S, P = load(exps)
    rows, out = [], [r"\begin{table*}[p]",
                     r"\caption{Control families: $\Omega$ of HyperCommon (HC), its triangles-only variant and the best "
                     r"other method under the oracle and label-free protocols, with the reference family each control "
                     r"changes (one assumption at a time; Sec.~\ref{sec:protocol}), mean over the noise range of each "
                     r"family ($p\in\{0,0.1,\dots,0.8\}$ for the controls, $p\in\{0,\dots,0.6\}$ for the swap controls "
                     r"and their references, $p\in[0,0.8]$ for the other references). The controls are compared with "
                     r"the methods that lead on the main benchmarks (DEMON, ANGEL, SCAN and SCAN-J are left out). "
                     r"Mark on HC as in Table~\ref{tab:labelfree}.}",
                     r"\label{tab:S-controls}", r"\scriptsize\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}",
                     r"\begin{tabular}{lcccccc|lcc}",
                     r" & \multicolumn{3}{c}{oracle} & \multicolumn{3}{c|}{label-free} & reference & HC oracle & HC free \\",
                     r"Control & HC & tri. & best other & HC & tri. & best other & & & \\", r"\hline"]
    for exp, fam, lab, rexp, rfam, rlab in CONTROLS:
        cs = [cell(S, P, exp, fam, prot) for prot in ("oracle", "label_free")]
        rs = [cell(S, P, rexp, rfam, prot) for prot in ("oracle", "label_free")]
        cells = []
        for c in cs:
            if c is None:
                cells += ["--"] * 3
                continue
            cells += [f"{c['hc']:.2f}$^{{{c['mark']}}}$", f"{c['tri']:.2f}", f"{c['best']:.2f} ({NAME[c['best_m']]})"]
        ref = [("--" if r is None else f"{r['hc']:.2f}") for r in rs]
        out.append(lab + " & " + " & ".join(cells) + f" & {rlab} & " + " & ".join(ref) + r" \\")
        for prot, c, r in zip(("oracle", "label_free"), cs, rs):
            if c:
                rows.append(dict(experiment=exp, family=fam, protocol=prot, **c,
                                 ref_hc=r["hc"] if r else np.nan, ref_best=r["best"] if r else np.nan,
                                 ref_best_m=r["best_m"] if r else ""))
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}", ""]
    # null graphs: what each method reports when there is nothing to find
    N = S[(S.experiment == "nulls") & (S.protocol == "label_free")]
    meths = [m for m in ["hc", "hc_tri", "linkcomm", "cpm", "leiden", "leiden_cpm", "infomap", "walktrap", "lpa",
                         "sbm", "slpa", "demon", "angel", "scan", "scan_j"] if m in set(N.method)] or ["hc"]
    out += [r"\begin{table*}[p]",
            r"\caption{Null graphs with triadic closure but no planted communities: number of communities and "
            r"coverage of each method's label-free cover (mean over ten graphs). RGG: random geometric graph on "
            r"the unit torus, 2000 nodes, mean degree 8 (clustering $0.59$); HK: Holme--Kim graph, 2000 nodes, "
            r"$m=4$, triangle probability $0.9$ (clustering $0.43$).}",
            r"\label{tab:S-nulls}", r"\scriptsize\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}",
            r"\begin{tabular}{ll" + "c" * len(meths) + "}",
            r"Graph & & " + " & ".join(NAME.get(m, m) for m in meths) + r" \\", r"\hline"]
    for fam, lab in (("rgg_null", "RGG"), ("hk_null", "HK")):
        for metric, mlab, fmt in (("n_comms", "communities", "%.0f"), ("coverage", "coverage", "%.2f")):
            g = N[(N.family == fam) & (N.metric == metric)].set_index("method")["mean"]
            out.append(f"{lab if metric == 'n_comms' else ''} & {mlab} & " +
                       " & ".join("--" if m not in g else fmt % g[m] for m in meths) + r" \\")
            for m in meths:
                if m in g:
                    rows.append(dict(experiment="nulls", family=fam, protocol="label_free", metric=metric, method=m,
                                     value=g[m]))
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    open(os.path.join(SEC, "supp_controls.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "controls3.csv"), index=False)
    print("controls written")


if __name__ == "__main__":
    main()
