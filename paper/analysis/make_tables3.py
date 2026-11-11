"""
Benchmark tables of the manuscript, generated from the result files so that
no number is typed by hand.

  sections/table_benchmarks.tex   label-free Omega, every ring setting (three overlap levels
                                  for every unequal family) and LFR family, all twelve other
                                  methods and the three ablations (tab:labelfree)
  sections/table_scorecard.tex    HC minus the best other method under the four protocols,
                                  with the paired-interval mark, and the counts (tab:scorecard)
  sections/supp_benchmarks.tex    Supplemental tables: oracle, frozen and native Omega; ONMI and
                                  all-community recovery of the oracle and label-free covers;
                                  small/large-community recovery; coverage and number of
                                  communities; the paired intervals behind the counts; the
                                  selected parameters; the grid-edge audit
  figures/data/scorecard3.csv     the counts and margins, for the text
Sources: bench/results2/<exp>_summary2.csv, _paired2.csv and _selected2.csv (analysis2.py).

    python make_tables3.py
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.join(HERE, "..", "bench", "results2")
SEC = os.path.join(HERE, "..", "sections")
DATA = os.path.join(HERE, "..", "figures", "data")

UNEQ = [("sizes_only_k8", r"equal $k$, unequal $s$, $k{=}8$"),
        ("dense_even", r"equal $k$, unequal $s$, $k{=}16$"),
        ("k_only_s100", r"equal $s$, unequal $k$"),
        ("two_worlds", "two worlds"),
        ("three_big_many_small", "three big, many small"),
        ("gradient10", "gradient")]
# (experiment, family, label, o_n, group)
ROWS = []
for k, s in ((8, 50), (16, 100)):
    for i, o in enumerate(("00", "05", "10")):
        ROWS.append(("rings_uniform", f"k{k}_s{s}_o{o}", rf"uniform, $k{{=}}{k}$" if i == 0 else "", o, f"u{k}"))
for fam, lab in UNEQ:
    for i, o in enumerate(("00", "05", "10")):
        exp, f = ("rings_hetero", fam) if o == "05" else ("rings_hetero_overlap", f"{fam}_o{o}")
        ROWS.append((exp, f, lab if i == 0 else "", o, fam))
for fam, lab in (("k8_s50_o05_preferential", r"uniform $k{=}8$, $\alpha{=}1$"),
                 ("k8_s50_o05_preferential2", r"uniform $k{=}8$, $\alpha{=}2$"),
                 ("two_worlds_preferential", r"two worlds, $\alpha{=}1$"),
                 ("two_worlds_preferential2", r"two worlds, $\alpha{=}2$")):
    ROWS.append(("rings_pref", fam, lab, "05", "pref"))
RING_ROWS = list(ROWS)
for exp, fam, lab in (("lfr", "on00", r"LFR, $o_n{=}0$"), ("lfr", "on10", r"LFR, $o_n{=}0.1$"),
                      ("lfr", "on30", r"LFR, $o_n{=}0.3$"), ("lfr_slices", "k10", r"LFR, $\langle k\rangle{=}10$"),
                      ("lfr_slices", "om4", r"LFR, $o_m{=}4$"), ("lfr_slices", "big", r"LFR, sizes 50--250"),
                      ("lfr_clustered", "on10_C03", r"LFR, $C_{\rm LFR}{=}0.3$"),
                      ("lfr_clustered", "on10_C05", r"LFR, $C_{\rm LFR}{=}0.5$")):
    ROWS.append((exp, fam, lab, "", "lfr"))

ABL = [("hc", "HC"), ("hc_econn", "e-c."), ("hc_tri", "tri."), ("hc_2test", "2-test")]
OTHERS = [("linkcomm", "Link"), ("cpm", "CPM"), ("leiden", "Leiden"), ("leiden_cpm", "L-CPM"),
          ("infomap", "Info."), ("walktrap", "Walk."), ("lpa", "LPA"), ("sbm", "SBM"), ("slpa", "SLPA"), ("oslom", "OSLOM"),
          ("demon", "DEMON"), ("angel", "ANGEL"), ("scan", "SCAN"), ("scan_j", "SCAN-J")]
METHODS = ABL + OTHERS
OTHER_IDS = [m for m, _ in OTHERS]
NAME = dict(METHODS)
PROTS = [("oracle", "oracle"), ("label_free", "label-free"), ("frozen", "frozen"), ("native", "native")]
HEAD2 = (r" & & \multicolumn{4}{c}{HC and ablations} & \multicolumn{2}{c}{relatives} &"
         r"\multicolumn{6}{c}{partitions} & \multicolumn{6}{c}{overlapping methods} \\")
LEGEND = (r"HC: HyperCommon; e-c.: its triples joined only through shared edges; tri.: triangles only "
          r"(J-thresholded 3-clique percolation~\cite{palla2005uncovering,kumpula2008sequential}); 2-test: "
          r"third pair not tested; Link: link communities with the same similarity; CPM: clique percolation; "
          r"Leiden: Leiden with the Reichardt--Bornholdt modularity; L-CPM: Leiden with the Constant Potts "
          r"Model; Info.: Infomap; Walk.: Walktrap, walks of length 4, dendrogram cut at maximum modularity; LPA: label propagation; "
          r"SBM: nested degree-corrected stochastic block model, each level of its hierarchy a parameter value; "
          r"OSLOM: order statistics local optimisation, significance threshold $t$ as its parameter; "
          r"SCAN: SCAN with cosine similarity; SCAN-J: SCAN with Jaccard similarity; $\alpha$: exponent of the "
          r"preferential (degree-biased) noise.")
# the SBM's frozen level can be missing on a test graph whose fitted hierarchy is shallower
FROZEN_DASH = (r" Dashes: the level of the stochastic block model's hierarchy frozen on the calibration "
               r"graphs does not exist on some test graphs of that family, whose fitted hierarchy has fewer "
               r"levels; on these rows the best other method is taken over the remaining methods, which cannot "
               r"change a mark, since the stochastic block model's oracle $\Omega$ there is far below that of "
               r"the best other method (Table~\ref{tab:oracle}).")
NOISE = (r"mean over the noise range ($p\in[0,0.8]$ for rings, $\mu\in[0.1,0.6]$ for LFR), ten test graphs "
         r"per noise level")
MARK = {"+": "^{+}", "=": "^{=}", "-": "^{-}"}


def load():
    exps = sorted({r[0] for r in ROWS})
    S = pd.concat([pd.read_csv(os.path.join(R2, f"{e}_summary2.csv")).assign(experiment=e) for e in exps],
                  ignore_index=True)
    P = pd.concat([pd.read_csv(os.path.join(R2, f"{e}_paired2.csv")).assign(experiment=e) for e in exps],
                  ignore_index=True)
    return S, P


def get(S, exp, fam, prot, metric):
    g = S[(S.experiment == exp) & (S.family == fam) & (S.protocol == prot) & (S.metric == metric)]
    return g.set_index("method")["mean"]


def verdict(S, P, exp, fam, prot, metric="omega"):
    """HC against the best and second-best other method (ablations excluded): margins with
    paired 95% intervals (HC minus other) and the mark +, = or -."""
    g = get(S, exp, fam, prot, metric)
    if g.empty or "hc" not in g:
        return None
    o = g.reindex(OTHER_IDS).dropna().sort_values(ascending=False)
    out = dict(experiment=exp, family=fam, protocol=prot, hc=g["hc"])
    for rank, m in enumerate(o.index[:2], 1):
        pr = P[(P.experiment == exp) & (P.family == fam) & (P.protocol == prot) & (P.metric == metric)
               & (P.method == m)]
        if pr.empty:
            continue
        pr = pr.iloc[0]
        out[f"m{rank}"], out[f"v{rank}"] = m, o[m]
        out[f"d{rank}"], out[f"lo{rank}"], out[f"hi{rank}"] = -pr["diff"], -pr.hi, -pr.lo
    out["mark"] = "+" if out["lo1"] > 0 else "-" if out["hi1"] < 0 else "="
    return out


def write(name, lines):
    open(os.path.join(SEC, name), "w", encoding="utf-8").write("\n".join(lines) + "\n")


def full_table(S, P, prot, metric, label, caption, rows=ROWS, marks=True, fmt="%.2f", star=True, bold=True):
    cols = "llcccc|cc|cccccc|" + "c" * 6        # overlapping: SLPA, OSLOM, DEMON, ANGEL, SCAN, SCAN-J
    if bold and "Bold:" not in caption:
        caption += r" Bold: best in row, within 0.005."
    out = [r"\begin{table*}" + ("[p]" if not star else "[t]"), r"\caption{" + caption + "}", r"\label{" + label + "}",
           r"\scriptsize\setlength{\tabcolsep}{1.35pt}",
           # twenty columns, scaled to the text width (the one-column Supplemental pages are narrower)
           r"\resizebox{\textwidth}{!}{%", r"\begin{tabular}{" + cols + "}", r"\hline\hline",
           HEAD2, r"Family & $o_n$ & " + " & ".join(n for _, n in METHODS) + r" \\", r"\hline"]
    last = None
    for exp, fam, lab, o, grp in rows:
        g = get(S, exp, fam, prot, metric)
        if g.empty:
            continue
        if last is not None and grp != last:
            out.append(r"\hline")
        last = grp
        vals = {m: g.get(m, np.nan) for m, _ in METHODS}
        best = np.nanmax(list(vals.values()))
        v = verdict(S, P, exp, fam, prot, metric) if marks else None
        cells = []
        for m, _ in METHODS:
            x = vals[m]
            if np.isnan(x):
                cells.append("--")
                continue
            s = fmt % x
            if m == "hc" and v is not None:
                s += "$" + MARK[v["mark"]] + "$"
            cells.append(r"\textbf{" + s + "}" if bold and abs(x - best) < 0.005 else s)
        on = f"{int(o) / 100:.2f}" if o else ""
        out.append(lab + " & " + on + " & " + " & ".join(cells) + r" \\")
    out += [r"\hline\hline", r"\end{tabular}}", r"\end{table*}"]
    return out


def counts(V, prot, rows):
    sub = [v for v in V if v and v["protocol"] == prot and (v["experiment"], v["family"]) in rows]
    c = {k: sum(v["mark"] == k for v in sub) for k in "+=-"}
    return c, len(sub)


def scorecard(S, P):
    V = [verdict(S, P, exp, fam, prot) for exp, fam, *_ in ROWS for prot, _ in PROTS]
    pd.DataFrame([v for v in V if v]).to_csv(os.path.join(DATA, "scorecard3.csv"), index=False)
    out = [r"\begin{table*}[t]",
           r"\caption{Scorecard of the benchmarks under the four protocols: $\Omega$ of HyperCommon minus that "
           r"of the best other method (ablations excluded), " + NOISE + r", with the best other method in "
           r"brackets when HyperCommon is not significantly ahead of it and its difference rounds to 0.00 or "
           r"below. Mark: the paired 95\% bootstrap interval of the difference over the ten test seeds "
           r"lies above zero ($+$), contains zero ($=$) or lies below zero ($-$); intervals and the second-best "
           r"method in Table~\ref{tab:S-intervals}. Oracle: best parameter per graph; label-free: largest "
           r"overlapping modularity, Eq.~\eqref{eq:qov}; frozen: one parameter per family from calibration "
           r"graphs; native: each method's own default or selection rule where it has one "
           r"(Sec.~\ref{sec:protocol}), the label-free choice otherwise. Counts are not corrected for "
           r"multiplicity. Method abbreviations as in Table~\ref{tab:labelfree}; the full tables are "
           r"Tables~\ref{tab:labelfree} and~\ref{tab:oracle}--\ref{tab:native}.}",
           r"\label{tab:scorecard}", r"\scriptsize\setlength{\tabcolsep}{3pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{llcccc}",
           r"Family & $o_n$ & " + " & ".join(n for _, n in PROTS) + r" \\", r"\hline"]
    idx = {(v["experiment"], v["family"], v["protocol"]): v for v in V if v}
    last = None
    for exp, fam, lab, o, grp in ROWS:
        if last is not None and grp != last:
            out.append(r"\hline")
        last = grp
        cells = []
        for prot, _ in PROTS:
            v = idx.get((exp, fam, prot))
            if v is None:
                cells.append("--")
                continue
            s = f"${v['d1']:+.2f}^{{{v['mark'] if v['mark'] != '-' else '-'}}}$"
            s = s.replace("^{+}", "^{+}").replace("^{=}", "^{=}")
            if v["d1"] < 0.005 and v["mark"] != "+":
                s += f" ({NAME[v['m1']]})"
            cells.append(s)
        on = f"{int(o) / 100:.2f}" if o else ""
        out.append(lab + " & " + on + " & " + " & ".join(cells) + r" \\")
    ring = {(e, f) for e, f, *_ in RING_ROWS}
    lfr = {(e, f) for e, f, *_ in ROWS if (e, f) not in ring}
    out.append(r"\hline")
    for nm, rs in (("rings", ring), ("LFR", lfr)):
        cells = []
        for prot, _ in PROTS:
            c, n = counts(V, prot, rs)
            cells.append(f"{c['+']}/{c['=']}/{c['-']} of {n}")
        out.append(rf"$+/{{=}}/-$, {nm} & & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    write("table_scorecard.tex", out)
    summ = []
    for nm, rs in (("rings", ring), ("lfr", lfr)):
        for prot, _ in PROTS:
            c, n = counts(V, prot, rs)
            summ.append(dict(set=nm, protocol=prot, ahead=c["+"], level=c["="], behind=c["-"], n=n))
    pd.DataFrame(summ).to_csv(os.path.join(DATA, "scorecard3_counts.csv"), index=False)
    print(pd.DataFrame(summ).to_string())
    return V


def intervals_table(V):
    out = [r"\begin{table*}[p]",
           r"\caption{The paired intervals behind Table~\ref{tab:scorecard}: $\Omega$ of HyperCommon minus that "
           r"of the best (first line of each cell) and second-best other method, with 95\% bootstrap intervals "
           r"over the ten test seeds, ablations excluded.}",
           r"\label{tab:S-intervals}", r"\tiny\setlength{\tabcolsep}{1pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{llcccc}", r"Family & $o_n$ & " + " & ".join(n for _, n in PROTS) + r" \\", r"\hline"]
    idx = {(v["experiment"], v["family"], v["protocol"]): v for v in V if v}
    for exp, fam, lab, o, grp in ROWS:
        cells = []
        for prot, _ in PROTS:
            v = idx.get((exp, fam, prot))
            if v is None or "m2" not in v:
                cells.append("--")
                continue
            # top-aligned, so the best line sits on the row's own baseline
            cells.append(r"\begin{tabular}[t]{@{}c@{}}" + f"{NAME[v['m1']]} ${v['d1']:+.3f}$ [${v['lo1']:+.3f}$, ${v['hi1']:+.3f}$]"
                         + r"\\" + f"{NAME[v['m2']]} ${v['d2']:+.3f}$ [${v['lo2']:+.3f}$, ${v['hi2']:+.3f}$]" + r"\end{tabular}")
        lab2 = lab if lab else r"\phantom{x}"
        on = f"{int(o) / 100:.2f}" if o else ""
        out.append(lab2 + " & " + on + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    return out


SEL_METHODS = [("hc", "HC"), ("linkcomm", "Link"), ("cpm", "CPM"), ("leiden", "Leiden"),
               ("leiden_cpm", "L-CPM"), ("infomap", "Infomap"), ("slpa", "SLPA")]


def selected():
    exps = sorted({r[0] for r in ROWS})
    return pd.concat([pd.read_csv(os.path.join(R2, f"{e}_selected2.csv"), dtype={"param": str}).assign(experiment=e)
                      for e in exps if os.path.exists(os.path.join(R2, f"{e}_selected2.csv"))], ignore_index=True)


def selector_table(S, D):
    out = [r"\begin{table*}[p]",
           r"\caption{The label-free selector against the oracle on the ring settings, for HyperCommon, link "
           r"communities, clique percolation, the two Leiden objectives, Infomap and SLPA: the ratio of "
           r"label-free to oracle $\Omega$ (both mean over $p\in[0,0.8]$) and the "
           r"median number of communities of the selected covers over the test graphs, label-free/oracle "
           r"(number of planted rings in the first column). For a partition the selector is the "
           r"Newman--Girvan modularity.}",
           r"\label{tab:S-selector}", r"\tiny\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{llc" + "cc" * len(SEL_METHODS) + "}",
           r" & & & " + " & ".join(r"\multicolumn{2}{c}{" + n + "}" for _, n in SEL_METHODS) + r" \\",
           r"Family & $o_n$ & rings & " + " & ".join(r"ratio & comm." for _ in SEL_METHODS) + r" \\", r"\hline"]
    import experiments as E
    for exp, fam, lab, o, grp in RING_ROWS:
        fdef = E.EXPERIMENTS[exp]["families"].get(fam)
        if fdef is None:
            continue
        cells = [str(len(fdef["sizes"]))]
        orc, lf = get(S, exp, fam, "oracle", "omega"), get(S, exp, fam, "label_free", "omega")
        if orc.empty:
            continue
        d = D[(D.experiment == exp) & (D.family == fam)]
        for m, _ in SEL_METHODS:
            if m not in orc or m not in lf:
                cells += ["--", "--"]
                continue
            nf = d[(d.method == m) & (d.protocol == "label_free")].n_comms.median()
            no = d[(d.method == m) & (d.protocol == "oracle")].n_comms.median()
            cells += [f"{lf[m] / orc[m]:.2f}", f"{nf:.0f}/{no:.0f}"]
        on = f"{int(o) / 100:.2f}" if o else ""
        out.append(lab + " & " + on + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    return out


def gridedge_table(D):
    import methods2
    rows = []
    for m, n in METHODS:
        grid = methods2.GRIDS2.get(m)
        if not grid or len(grid) < 2 or grid[0] is None:
            continue
        if isinstance(grid[0], tuple):
            es = sorted({g[0] for g in grid})
            first = lambda p: float(str(p).strip("()").split(",")[0])
            low = lambda p, es=es: first(p) == es[0]
            high = lambda p, es=es: first(p) == es[-1]
        else:
            lo, hi = min(grid), max(grid)
            low = lambda p, lo=lo: float(p) == lo
            high = lambda p, hi=hi: float(p) == hi
        cells = []
        for prot, _ in PROTS[:3]:
            d = D[(D.method == m) & (D.protocol == prot)]
            if not len(d):
                cells.append("--")
                continue
            a, b = d.param.map(low).mean(), d.param.map(high).mean()
            cells.append(f"{100 * (a + b):.0f}\\% ({100 * a:.0f}/{100 * b:.0f})")
        rows.append(n + " & " + " & ".join(cells) + r" \\")
    return [r"\begin{table}[p]",
            r"\caption{Grid-edge audit: share of the test graphs (all ring and LFR families) whose selected "
            r"parameter lies at an end of the method's grid (Table~\ref{tab:grids}); for SCAN the end points of "
            r"$\varepsilon$; in brackets the shares at the low and the high end. HyperCommon's lowest "
            r"threshold, $0.02$, is an end of its grid, and selecting it means the covers coalesce. The low end "
            r"of clique percolation, clique size 3, is the smallest possible.}",
            r"\label{tab:S-gridedge}", r"\footnotesize", r"\begin{ruledtabular}", r"\begin{tabular}{lccc}",
            r"Method & oracle & label-free & frozen \\", r"\hline"] + rows + \
           [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]


def main():
    import sys
    sys.path.insert(0, os.path.join(HERE, "..", "bench"))
    S, P = load()
    V = scorecard(S, P)

    def ci(prot):
        om = S[(S.metric == "omega") & (S.protocol == prot)]
        hw = ((om.hi - om.lo) / 2).dropna()
        return (rf"95\% bootstrap half-widths of these means over the ten test seeds have median "
                rf"{hw.median():.3f} and maximum {hw.max():.3f}.")

    lf = full_table(S, P, "label_free", "omega", "tab:labelfree",
                    r"Accuracy on the benchmarks as a user would run the methods, every parameter chosen without "
                    r"labels by one rule: on each test graph, the parameter whose cover has the largest overlapping "
                    r"modularity, Eq.~\eqref{eq:qov}, which for a partition is the Newman--Girvan modularity "
                    r"(label-free protocol). $\Omega$, " + NOISE + ". " + LEGEND + r" Parameter-free methods "
                    r"(Walktrap, LPA) are as under every protocol. Superscripts on HC: the paired 95\% bootstrap "
                    r"interval of its difference from the best other method (ablations excluded) lies above zero "
                    r"($+$), contains zero ($=$) or lies below zero ($-$). Bold: best in row, within 0.005. "
                    + ci("label_free"))
    write("table_benchmarks.tex", lf)
    sup = []
    for prot, label, txt in (
            ("oracle", "tab:oracle", r"every method at its best parameter on each test graph (oracle protocol)"),
            ("frozen", "tab:frozen", r"each method's parameter frozen on three calibration graphs per noise level "
                                     r"of the same family, all noise levels pooled (frozen protocol)"),
            ("native", "tab:native", r"each method at its own default or selection rule: Leiden at "
                                     r"$\gamma=1$, L-CPM at $\gamma$ equal to the edge density, link communities "
                                     r"at maximum partition density, Infomap at Markov time 1, SLPA at $r=0.1$, OSLOM at $t=0.1$, "
                                     r"SBM at the finest level of its hierarchy, and the label-free choice for the methods without one (native protocol)")):
        sup += full_table(S, P, prot, "omega", label, r"As Table~\ref{tab:labelfree}, with " + txt + ". " + ci(prot)
                          + (FROZEN_DASH if prot == "frozen" else ""),
                          star=False) + [""]
    for prot, nm in (("oracle", "oracle"), ("label_free", "label-free")):
        sup += full_table(S, P, prot, "onmi", f"tab:onmi-{nm}",
                          rf"ONMI (McDaid \emph{{et al.}}, max normalisation) of the {nm} covers of "
                          rf"Table~\ref{{{'tab:oracle' if prot == 'oracle' else 'tab:labelfree'}}}; superscripts on "
                          r"HC as in Table~\ref{tab:labelfree}, for ONMI.", star=False) + [""]
        sup += full_table(S, P, prot, "rec_mean", f"tab:S-allcomm-{nm}",
                          rf"Per-community recovery over all true communities (mean over true communities of the "
                          rf"Jaccard similarity with the best-matching detected community) of the {nm} covers; "
                          r"superscripts as in Table~\ref{tab:labelfree}, for this score.", star=False) + [""]
    uneq = [r for r in RING_ROWS if r[4] in dict(UNEQ)]
    for metric, label, txt in (("rec_small", "tab:recsmall", "no larger than"),
                               ("rec_large", "tab:reclarge", "larger than")):
        sup += full_table(S, P, "oracle", metric, label,
                          rf"Per-community recovery of the true communities {txt} the median size, oracle covers, "
                          r"on the unequal-ring families (equal $s$ with unequal $k$ has no community larger than "
                          r"the median). " + (r"Numbers of rings no larger than the median size: $20$, $20$, $20$, "
                          r"$20$, $16$, $5$ (families in order)." if metric == "rec_small" else
                          r"Numbers of rings larger than the median size: $10$, $10$, $5$, $3$, $5$ (families "
                          r"in order)."), rows=uneq, marks=False,
                          star=False) + [""]
    for prot, nm in (("oracle", "oracle"), ("label_free", "label-free")):
        sup += full_table(S, P, prot, "coverage", f"tab:S-coverage-{nm}",
                          rf"Coverage (fraction of nodes in at least one community) of the {nm} covers.",
                          marks=False, star=False, bold=False) + [""]
        sup += full_table(S, P, prot, "n_comms", f"tab:S-ncomms-{nm}",
                          rf"Number of communities of the {nm} covers (planted: Table~\ref{{tab:shapes}}).",
                          marks=False, fmt="%.0f", star=False, bold=False) + [""]
    sup += intervals_table(V) + [""]
    D = selected()
    sup += selector_table(S, D) + [""] + gridedge_table(D)
    write("supp_benchmarks.tex", sup)


if __name__ == "__main__":
    main()
