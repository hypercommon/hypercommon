"""
Real-network tables of the manuscript.

  sections/table_real.tex        SNAP Amazon, DBLP, YouTube (full graphs): F1 and Omega at
                                 the best and the label-free parameter, coverage at the
                                 label-free parameter, and Leiden and Infomap at their
                                 defaults (tab:real)
  sections/table_real_more.tex   fourteen smaller networks, F1 and Omega (tab:realmore)
  sections/supp_real.tex         Supplemental: every method on the fourteen networks
                                 (tab:S-realnet-f1, tab:S-realnet-omega) and F1 of the SNAP
                                 covers against the distinct ground-truth sets
Sources: bench/results/real_big_scores.csv (score_big.py), bench/results/real_more_v3.csv.

    python make_tables3_real.py
"""
import os

import numpy as np
import pandas as pd

import make_tables2 as T2

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "bench", "results")
SEC = os.path.join(HERE, "..", "sections")

REAL_METHODS = [("components", "components"), ("hc", "HC"), ("hc_family", r"HC$_f$"),
                ("hc_econn", "e-conn."), ("hc_tri", "tri."), ("hc_2test", "2-test"),
                ("linkcomm", "Link"), ("leiden", "Leiden"), ("leiden@1.0", r"Leiden, $\gamma{=}1$"),
                ("leiden_cpm", "L-CPM"), ("infomap", "Infomap"), ("infomap@1.0", "Infomap, Markov time 1"),
                ("lpa", "LPA"), ("oslom", "OSLOM")]
NETS = [("amazon", "Amazon"), ("dblp", "DBLP"), ("youtube", "YouTube")]


def write(name, lines):
    open(os.path.join(SEC, name), "w", encoding="utf-8").write("\n".join(lines) + "\n")


def picks(setting, f1col="f1"):
    d = pd.read_csv(os.path.join(RES, "real_big_scores.csv"), dtype={"param": str})
    d = d[d.setting == setting].copy()
    fam = d[d.method.str.startswith("hcf_")].copy()
    fam["method"] = "hc_family"
    fixed = [d[(d.method == m) & (d.param == "1.0")].assign(method=f"{m}@1.0") for m in ("leiden", "infomap")]
    d = pd.concat([d[~d.method.str.startswith("hcf_")], fam] + fixed, ignore_index=True)
    out = {}
    for (net, m), g in d.groupby(["network", "method"]):
        out[(net, m)] = dict(best_f1=g.loc[g[f1col].idxmax()], best_om=g.loc[g.omega_dedup.idxmax()],
                             free=g.loc[g["eq"].idxmax()])
    return out


def real_table(setting="full", f1col="f1", label="tab:real", caption=None):
    P = picks(setting, f1col)
    cols = [("best_f1", f1col), ("free", f1col), ("best_om", "omega_dedup"), ("free", "omega_dedup"),
            ("free", "coverage")]
    rows = []
    for m, lab in REAL_METHODS:
        vals = []
        for net, _ in NETS:
            p = P.get((net, m))
            vals += [float(p[pk][c]) if p is not None else np.nan for pk, c in cols]
        if not np.all(np.isnan(vals)):
            rows.append((m, lab, vals))
    ncol = len(cols) * len(NETS)
    best = [np.nanmax([v[j] for m, _, v in rows if m != "components"] or [np.nan]) for j in range(ncol)]
    out = [r"\begin{table*}[t]", r"\caption{" + caption + "}", r"\label{" + label + "}",
           r"\scriptsize\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{l" + "c" * ncol + "}",
           " & " + " & ".join(r"\multicolumn{5}{c}{" + n + "}" for _, n in NETS) + r" \\",
           " & " + " & ".join([r"F1 best & F1 free & $\Omega$ best & $\Omega$ free & cov."] * 3) + r" \\",
           r"\hline"]
    for m, lab, vals in rows:
        cc = []
        for j, v in enumerate(vals):
            if np.isnan(v):
                cc.append("--")
                continue
            s = ("$%.3f$" % v).replace("$-", "$-") if v < 0 else "%.3f" % v
            if j % 5 == 4:
                s = "%.2f" % v
            elif m != "components" and abs(v - best[j]) < 5e-4:
                s = r"\textbf{" + s + "}"
            cc.append(s)
        out.append(lab + " & " + " & ".join(cc) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    return out


CAP_FULL = (r"Detection on the full SNAP Amazon (334\,863 nodes), DBLP (317\,080) and YouTube (1\,134\,890) "
            r"graphs, scored on the nodes of their ground-truth communities (top 5000 for Amazon and DBLP, all "
            r"groups of at least three nodes for YouTube), to which every detected community is restricted. "
            r"Best-match F1 against the recorded communities and $\Omega$ against their distinct node sets, at "
            r"each method's best parameter by that metric (best) and at the parameter chosen without labels by "
            r"maximal overlapping modularity on the full graph (free); cov.: fraction of the scored nodes in at "
            r"least one community at the free parameter. HC$_f$: the seven commonality functions together, the "
            r"free choice taken over all of them. Leiden at $\gamma=1$ and Infomap at Markov time 1 are the two "
            r"methods at their defaults (one parameter, so best and free coincide). Connected components are a "
            r"trivial reference. On YouTube HyperCommon uses degree-ratio pruning, exact for every $t\ge0.02$ "
            r"(Sec.~\ref{sec:cost}), which excludes link communities and the other functions; clique percolation "
            r"and SLPA are compared on the induced subgraphs (Table~\ref{tab:realinduced}), and SCAN, Walktrap, "
            r"DEMON and ANGEL were not run on these graphs, to bound the computing time of the large runs. "
            r"Stochastic methods: one seed. F1 against the distinct node sets is in "
            r"Table~\ref{tab:S-realnet-dedup}; it changes F1 by at most 0.03 (at most 0.01 for HyperCommon and "
            r"its variants) and reverses the order of no two methods in a column. The top-5000 lists of Amazon and DBLP "
            r"frequently consist of near-cliques~\cite{chen2024parallel}, which can favour methods that find "
            r"small dense groups. Bold: best in column.")


def supp_real():
    out = real_table("full", "f1_dedup", "tab:S-realnet-dedup",
                     r"As Table~\ref{tab:real}, with F1 scored against the distinct node sets of the recorded "
                     r"communities rather than against the top-5000 lists, which repeat some sets.")
    raw = pd.read_csv(os.path.join(RES, "real_more_v3.csv"), dtype={"param": str})
    d = raw[raw.error.isna()].copy()
    key = ["network", "method", "seed"]
    meths = ["hc", "hc_econn", "hc_tri", "hc_2test"] + T2.MORE_OTHERS
    names = dict(T2.METHODS) | {"walktrap": "Walktrap", "scan_j": "SCAN-J", "oslom": "OSLOM"}
    for metric, mname in (("f1", "best-match F1"), ("omega", r"$\Omega$"), ("coverage", "coverage")):
        for prot in ("best", "free"):
            col = metric if prot == "best" else "eq"
            if metric == "coverage" and prot == "best":
                continue
            sel = d.loc[d.groupby(key)[col if metric != "coverage" else "eq"].idxmax()] if prot == "free" else \
                d.loc[d.groupby(key)[metric].idxmax()]
            T = sel.groupby(["network", "method"])[metric].mean().unstack()
            T = T.reindex(columns=[m for m in meths if m in T.columns])
            cap = (rf"Every method on the fourteen smaller networks: {mname} at the "
                   + ("best parameter" if prot == "best" else "label-free parameter")
                   + r" (stochastic methods: mean over five seeds)"
                   + (r"." if metric == "coverage" else
                      r", with each method's rank among the methods shown (ablations included) in the last row "
                      r"as the mean rank over the fourteen networks. Bold: best in row, within 0.005."))
            lab = f"tab:S-realnet-{metric}-{prot}"
            head = "Network & " + " & ".join(names.get(m, m) for m in T.columns) + r" \\"
            body = []
            ranks = T.rank(axis=1, ascending=False)
            for net, nlab in T2.MORE_NETS:
                if net not in T.index:
                    continue
                r = T.loc[net]
                top = r.max()
                cells = []
                for m in T.columns:
                    v = r[m]
                    # true minus sign, and no "-0.00"
                    s = "--" if np.isnan(v) else ("$-$%.2f" % -v if round(v, 2) < 0 else "%.2f" % abs(v)
                                                  if round(v, 2) == 0 else "%.2f" % v)
                    if metric != "coverage" and not np.isnan(v) and abs(v - top) < 5e-3:
                        s = r"\textbf{" + s + "}"
                    cells.append(s)
                body.append(nlab + " & " + " & ".join(cells) + r" \\")
            if metric != "coverage":
                body += [r"\hline", "mean rank & " + " & ".join("%.1f" % ranks[m].mean() for m in T.columns) + r" \\"]
            out += ["", r"\begin{table*}[p]", r"\caption{" + cap + "}", r"\label{" + lab + "}",
                    r"\tiny\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}",
                    r"\begin{tabular}{l" + "c" * len(T.columns) + "}", head, r"\hline"] + body + \
                   [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    write("supp_real.tex", out)


if __name__ == "__main__":
    write("table_real.tex", real_table("full", "f1", "tab:real", CAP_FULL))
    T2.real_table("induced")
    T2.more_table()
    supp_real()
    print("real tables written")
