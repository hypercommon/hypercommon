"""
Tables of the second benchmark round, generated from the result files so that no
number is typed by hand.

  sections/table_benchmarks.tex   mean Omega over the noise range: oracle (tab:oracle)
                                  and label-free (tab:labelfree) protocols
  sections/table_frozen.tex       the same, frozen protocol (appendix)
  sections/table_onmi.tex         ONMI of the oracle covers (appendix)
  sections/table_recovery.tex     per-community recovery of small / large communities (appendix)
  sections/table_overlap.tex      unequal rings at overlap 0 and 0.10 (appendix)
  sections/table_robust.tex       robustness to noise on the ring families
  sections/table_real.tex         SNAP Amazon, DBLP, YouTube, full graphs
  sections/table_real_induced.tex the same on the subgraphs induced by the labelled nodes
  sections/table_real_more.tex    fourteen smaller networks, F1 and Omega
  sections/table_family.tex       the commonality family on the benchmarks (family2)
  sections/table_geometry.tex     internal geometry sweep (rings_geometry)
Sources: bench/results2/<exp>_summary2.csv (analysis2.py), results2/robust_p.csv
(robust_p.py), bench/results/real_big_scores.csv (score_big.py), real_more_v3.csv.

    python make_tables2.py [benchmarks|robust|real|more|family|geometry ...]
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "..", "bench")
R2 = os.path.join(BENCH, "results2")
RES = os.path.join(BENCH, "results")
SEC = os.path.join(HERE, "..", "sections")
DATA = os.path.join(HERE, "..", "figures", "data")

ROWS = [  # (experiment, family, label)
    ("rings_uniform", "k8_s50_o00", r"uniform, $k{=}8$, $o_n{=}0$"),
    ("rings_uniform", "k8_s50_o05", r"uniform, $k{=}8$, $o_n{=}0.05$"),
    ("rings_uniform", "k8_s50_o10", r"uniform, $k{=}8$, $o_n{=}0.10$"),
    ("rings_uniform", "k16_s100_o00", r"uniform, $k{=}16$, $o_n{=}0$"),
    ("rings_uniform", "k16_s100_o05", r"uniform, $k{=}16$, $o_n{=}0.05$"),
    ("rings_uniform", "k16_s100_o10", r"uniform, $k{=}16$, $o_n{=}0.10$"),
    ("rings_hetero", "sizes_only_k8", r"equal $k$, unequal $s$ ($k{=}8$)"),
    ("rings_hetero", "dense_even", r"equal $k$, unequal $s$ ($k{=}16$)"),
    ("rings_hetero", "k_only_s100", r"equal $s$, unequal $k$"),
    ("rings_hetero", "two_worlds", "two worlds"),
    ("rings_hetero", "three_big_many_small", "three big, many small"),
    ("rings_hetero", "gradient10", "gradient of ten rings"),
    ("rings_pref", "k8_s50_o05_preferential", r"uniform $k{=}8$, pref.\ $\alpha{=}1$"),
    ("rings_pref", "k8_s50_o05_preferential2", r"uniform $k{=}8$, pref.\ $\alpha{=}2$"),
    ("rings_pref", "two_worlds_preferential", r"two worlds, pref.\ $\alpha{=}1$"),
    ("rings_pref", "two_worlds_preferential2", r"two worlds, pref.\ $\alpha{=}2$"),
    ("lfr", "on00", r"LFR, $o_n{=}0$"),
    ("lfr", "on10", r"LFR, $o_n{=}0.1$"),
    ("lfr", "on30", r"LFR, $o_n{=}0.3$"),
    ("lfr_slices", "k10", r"LFR, $\langle k\rangle{=}10$"),
    ("lfr_slices", "om4", r"LFR, $o_m{=}4$"),
    ("lfr_slices", "big", r"LFR, sizes 50--250"),
    ("lfr_clustered", "on10_C03", r"LFR, $C{=}0.3$"),
    ("lfr_clustered", "on10_C05", r"LFR, $C{=}0.5$"),
]
OVERLAP_ROWS = [("rings_hetero_overlap", f"{f}_o{o}", lab if o == "00" else "", o)
                for f, lab in (("two_worlds", "two worlds"), ("three_big_many_small", "three big, many small"),
                               ("gradient10", "gradient of ten rings"),
                               ("sizes_only_k8", r"equal $k$, unequal $s$ ($k{=}8$)"))
                for o in ("00", "10")]
METHODS = [("hc", "HC"), ("hc_econn", "e-conn."), ("hc_tri", "tri."), ("hc_2test", "2-test"),
           ("linkcomm", "Link"), ("cpm", "CPM"), ("leiden", "Leiden"), ("leiden_cpm", "L-CPM"),
           ("infomap", "Infomap"), ("lpa", "LPA"), ("slpa", "SLPA"), ("angel", "ANGEL"),
           ("demon", "DEMON"), ("scan", "SCAN")]
HC_VARIANTS = ("hc", "hc_econn", "hc_tri", "hc_2test")
NAME = dict(METHODS) | {"walktrap": "Walktrap", "scan_j": "SCAN-J", "components": "comp."}
NOISE = (r"mean over the noise range ($p\in[0,0.8]$ for rings, $\mu\in[0.1,0.6]$ for LFR), "
         r"ten test graphs per noise level")
LEGEND = (r"HC: HyperCommon; e-conn.: its triples joined only through shared edges; tri.: triangles only "
          r"(Jaccard-weighted 3-clique percolation); 2-test: third pair not tested; Link: link communities "
          r"with the same similarity; L-CPM: Leiden with the Constant Potts Model.")


def write(name, lines):
    open(os.path.join(SEC, name), "w", encoding="utf-8").write("\n".join(lines) + "\n")


def summaries():
    exps = sorted({r[0] for r in ROWS + OVERLAP_ROWS})
    return pd.concat([pd.read_csv(os.path.join(R2, f"{e}_summary2.csv")) for e in exps], ignore_index=True)


def bench_table(S, protocol, metric, label, caption, rows=ROWS, extra_col=None):
    sub = S[(S.protocol == protocol) & (S.metric == metric)]
    cols = "l" + ("c" if extra_col else "") + "c" * len(METHODS)
    head = "Family & " + (extra_col + " & " if extra_col else "") + " & ".join(m[1] for m in METHODS) + r" \\"
    out = [r"\begin{table*}[t]", r"\caption{" + caption + "}", r"\label{" + label + "}",
           r"\scriptsize\setlength{\tabcolsep}{2pt}", r"\begin{ruledtabular}", r"\begin{tabular}{" + cols + "}",
           head, r"\hline"]
    last = None
    for row in rows:
        exp, fam, lab = row[:3]
        g = sub[(sub.experiment == exp) & (sub.family == fam)].set_index("method")["mean"]
        if g.empty:
            continue
        grp = "lfr" if exp.startswith("lfr") else exp
        if last is not None and grp != last:
            out.append(r"\hline")
        last = grp
        vals = {m: g.get(m, np.nan) for m, _ in METHODS}
        best = np.nanmax(list(vals.values()))
        cells = ["--" if np.isnan(v) else (r"\textbf{%.2f}" % v if abs(v - best) < 0.005 else "%.2f" % v)
                 for v in vals.values()]
        lead = lab + (" & " + f"{int(row[3]) / 100:.2f}" if extra_col else "")
        out.append(lead + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    return out


def benchmarks():
    S = summaries()
    S.to_csv(os.path.join(DATA, "benchmarks_summary2.csv"), index=False)
    om = S[(S.metric == "omega")]
    hw = ((om.hi - om.lo) / 2).dropna()
    ci = (rf"95\% bootstrap half-widths over the test seeds have median {hw.median():.3f} and maximum "
          rf"{hw.max():.3f} (all intervals in the supplementary data).")
    t1 = bench_table(S, "oracle", "omega", "tab:oracle",
                     r"Accuracy on the benchmarks with every method at its best parameter on each test graph "
                     r"(oracle protocol): $\Omega$, " + NOISE + ". " + LEGEND +
                     r" Bold: best in row, within 0.005. " + ci)
    t2 = bench_table(S, "label_free", "omega", "tab:labelfree",
                     r"As Table~\ref{tab:oracle}, with every method's parameter chosen without labels by one "
                     r"rule: on each test graph, the parameter whose cover has the largest overlapping "
                     r"modularity (label-free protocol). Parameter-free methods (LPA) are as in "
                     r"Table~\ref{tab:oracle}.")
    write("table_benchmarks.tex", t1 + [""] + t2)
    write("table_frozen.tex", bench_table(
        S, "frozen", "omega", "tab:frozen",
        r"As Table~\ref{tab:oracle}, with each method's parameter frozen on three calibration graphs per noise "
        r"level of the same family, all noise levels pooled, and applied to the test graphs (frozen protocol)."))
    write("table_onmi.tex", bench_table(
        S, "oracle", "onmi", "tab:onmi",
        r"As Table~\ref{tab:oracle}, scored with ONMI (McDaid \emph{et al.}, max normalisation) on the same "
        r"covers, i.e.\ at each method's $\Omega$-optimal parameter."))
    hetero = [r for r in ROWS if r[0] == "rings_hetero"]
    write("table_recovery.tex",
          bench_table(S, "oracle", "rec_small", "tab:recsmall",
                      r"Per-community recovery on the unequal-ring families: mean over true communities no "
                      r"larger than the median size of the Jaccard similarity with the best-matching detected "
                      r"community, at each method's $\Omega$-optimal parameter (covers of "
                      r"Table~\ref{tab:oracle}), " + NOISE + ".", rows=hetero) + [""] +
          bench_table(S, "oracle", "rec_large", "tab:reclarge",
                      r"As Table~\ref{tab:recsmall}, for the true communities larger than the median size.",
                      rows=hetero))
    write("table_overlap.tex", bench_table(
        S, "oracle", "omega", "tab:overlaplevels",
        r"Unequal rings at overlap $o_n=0$ and $0.10$ (Table~\ref{tab:oracle} uses $0.05$): $\Omega$ under the "
        r"oracle protocol, " + NOISE + ".", rows=OVERLAP_ROWS, extra_col="$o_n$") + [""] + bench_table(
        S, "label_free", "omega", "tab:overlaplevelsfree",
        r"As Table~\ref{tab:overlaplevels}, label-free protocol.", rows=OVERLAP_ROWS, extra_col="$o_n$"))
    return S


def robust():
    d = pd.read_csv(os.path.join(R2, "robust_p.csv"))
    out = [r"\begin{table*}[t]",
           r"\caption{Robustness to noise on the ring families: mean $\Omega$ over the high-noise half "
           r"($p\in[0.4,0.8]$) and $p_{1/2}$, the rewired fraction at which the mean $\Omega$ over the ten test "
           r"graphs first falls below $0.5$ (linear interpolation; $>0.8$: never), for HyperCommon (HC) and the "
           r"best of the other methods by each measure, under the oracle and label-free protocols. Bold: HC "
           r"ahead or within 0.005 (0.01 for $p_{1/2}$).}",
           r"\label{tab:robust}", r"\footnotesize\setlength{\tabcolsep}{3pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{lcccccccc}",
           r" & \multicolumn{4}{c}{oracle} & \multicolumn{4}{c}{label-free} \\",
           r"Family & HC $\bar\Omega_{\rm hi}$ & best other & HC $p_{1/2}$ & best other & "
           r"HC $\bar\Omega_{\rm hi}$ & best other & HC $p_{1/2}$ & best other \\", r"\hline"]
    fmt_p = lambda v: "$>0.8$" if not np.isfinite(v) else f"{v:.2f}"
    last = None
    for exp, fam, lab in ROWS:
        if not exp.startswith("rings"):
            continue
        if last is not None and exp != last:
            out.append(r"\hline")
        last = exp
        cells = []
        for prot in ("oracle", "label_free"):
            g = d[(d.experiment == exp) & (d.family == fam) & (d.protocol == prot)].set_index("method")
            o = g.drop([m for m in HC_VARIANTS if m in g.index])
            for col, tol, f in (("high_noise_mean", 0.005, "{:.2f}".format), ("p_half", 0.01, fmt_p)):
                hv, bm = g.loc["hc", col], o[col].idxmax()
                bv = o.loc[bm, col]
                hs = f(hv)
                cells += [r"\textbf{" + hs + "}" if hv >= bv - tol else hs, f"{f(bv)} ({NAME.get(bm, bm)})"]
        out.append(lab + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    write("table_robust.tex", out)


GEOM_BASE = [("k8_s50_o05", r"uniform, $k{=}8$"), ("two_worlds", "two worlds"),
             ("three_big_many_small", "three big, many small"), ("gradient10", "gradient of ten rings")]
GEOM_Q = ("000", "010", "030", "060", "100")


def geometry():
    """Internal geometry sweep (rings_geometry): a fraction q of each ring's edges rewired
    inside the ring before the noise, from lattice rings (q=0) to random graphs (q=1)."""
    S = pd.read_csv(os.path.join(R2, "rings_geometry_summary2.csv"))
    P = pd.read_csv(os.path.join(R2, "rings_geometry_paired2.csv"))
    S2 = S[S.metric == "coverage"]
    S, P = S[S.metric == "omega"], P[P.metric == "omega"]
    out = [r"\begin{table*}[t]",
           r"\caption{Internal geometry: before the noise, a fraction $q$ of each ring's edges is rewired "
           r"inside the ring, from lattice rings ($q=0$) to rings that are random graphs ($q=1$), at fixed "
           r"size and mean degree; $o_n=0.05$. Uniform: 40 rings of 50 nodes, $k=8$; two worlds: five rings "
           r"of 200 nodes ($k=32$) and twenty of 50 ($k=8$); three big, many small: three of 400 ($k=32$) "
           r"and sixteen of 50 ($k=8$); gradient: ten rings of 600 to 30 nodes ($k=16$ to 8). "
           r"Mean $\Omega$ over $p\in\{0,0.1,\dots,0.6\}$, ten test graphs per noise "
           r"level, under the oracle and label-free protocols, for HyperCommon (HC), its triangles-only "
           r"variant and the best other method; $C$: mean clustering of the noiseless graphs (ten seeds); "
           r"cov.: coverage of HC's label-free cover. Superscripts on HC: paired 95\% interval of its "
           r"difference from the best other method excludes zero ($^+$ HC ahead, $^-$ behind). Bold: HC best or "
           r"within 0.005 of the best other method. Means are over $p\in\{0,0.1,\dots,0.6\}$, so the $q=0$ "
           r"entries differ from those of Tables~\ref{tab:labelfree} and~\ref{tab:oracle} for the same "
           r"families, which average over "
           r"$p\in[0,0.8]$.}",
           r"\label{tab:geometry}", r"\footnotesize\setlength{\tabcolsep}{3pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{lcccccccccc}",
           r" & & & \multicolumn{3}{c}{oracle} & \multicolumn{4}{c}{label-free} \\",
           r"Family & $q$ & $C$ & HC & tri. & best other & HC & tri. & best other & cov. \\", r"\hline"]
    import networkx as nx
    sys.path.insert(0, BENCH)
    import experiments as E
    X = E.EXPERIMENTS["rings_geometry"]
    rows = []
    for i, (base, lab) in enumerate(GEOM_BASE):
        if i:
            out.append(r"\hline")
        for j, q in enumerate(GEOM_Q):
            fam = f"{base}_q{q}"
            cells = []
            for prot in ("oracle", "label_free"):
                g = S[(S.family == fam) & (S.protocol == prot)].set_index("method")["mean"]
                if g.empty or "hc" not in g or "hc_tri" not in g:   # run not finished
                    cells += ["--"] * 3
                    continue
                o = g.drop([m for m in HC_VARIANTS if m in g.index])
                bm = o.idxmax()
                pr = P[(P.family == fam) & (P.protocol == prot) & (P.method == bm)].iloc[0]
                sig = "-" if pr.lo > 0 else "+" if pr.hi < 0 else ""
                hs = f"{g['hc']:.2f}" + (f"$^{sig}$" if sig else "")
                cells += [r"\textbf{" + hs + "}" if g["hc"] >= o[bm] - 0.005 else hs, f"{g['hc_tri']:.2f}",
                          f"{o[bm]:.2f} ({NAME.get(bm, bm)})"]
                rows.append(dict(family=base, q=int(q) / 100, protocol=prot, hc=g["hc"], hc_tri=g["hc_tri"],
                                 best_other=bm, best_other_mean=o[bm], diff=pr["diff"], lo=pr.lo, hi=pr.hi))
            clus = np.mean([nx.average_clustering(X["build"](X["families"][fam], 0.0, sd)[0]) for sd in range(10)])
            cv = S2[(S2.family == fam) & (S2.protocol == "label_free") & (S2.method == "hc")]["mean"]
            cov = f"{cv.iloc[0]:.2f}" if len(cv) else "--"
            cells.append(cov)
            out.append((lab if j == 0 else "") + f" & {int(q) / 100:.1f} & {clus:.2f} & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    write("table_geometry.tex", out)
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "geometry_summary2.csv"), index=False)


# ---------------------------------------------------------------------------- real networks
REAL_METHODS = [("components", "components"), ("hc", "HyperCommon"), ("hc_family", r"HC$_f$"),
                ("hc_econn", "edge-connected"), ("hc_tri", "triangles only"), ("hc_2test", "two-test"),
                ("linkcomm", "link comm."), ("cpm", "CPM"), ("leiden", "Leiden"), ("leiden_cpm", "L-CPM"),
                ("infomap", "Infomap"), ("lpa", "label prop."), ("slpa", "SLPA")]
REAL_NETS = [("amazon", "Amazon"), ("dblp", "DBLP"), ("youtube", "YouTube")]


def real_picks(setting):
    d = pd.read_csv(os.path.join(RES, "real_big_scores.csv"), dtype={"param": str})
    d = d[d.setting == setting].copy()
    fam = d[d.method.str.startswith("hcf_")].copy()
    fam["method"] = "hc_family"
    d = pd.concat([d[~d.method.str.startswith("hcf_")], fam], ignore_index=True)
    picks = {}
    for (net, m), g in d.groupby(["network", "method"]):
        picks[(net, m)] = dict(best_f1=g.loc[g.f1.idxmax()], best_om=g.loc[g.omega_dedup.idxmax()],
                               free=g.loc[g["eq"].idxmax()])
    return picks


def real_table(setting):
    P = real_picks(setting)
    cols = [("best_f1", "f1"), ("free", "f1"), ("best_om", "omega_dedup"), ("free", "omega_dedup")]
    rows = []
    for m, lab in REAL_METHODS:
        cells, vals = [], []
        for net, _ in REAL_NETS:
            p = P.get((net, m))
            for pick, col in cols:
                v = float(p[pick][col]) if p is not None else np.nan
                vals.append(v)
        if not np.all(np.isnan(vals)):
            rows.append((m, lab, vals))
    best = [np.nanmax([v[j] for m, _, v in rows if m != "components"] or [np.nan]) for j in range(len(rows[0][2]))]
    if setting == "full":
        cap = (r"Detection on the full SNAP Amazon (334\,863 nodes), DBLP (317\,080) and YouTube (1\,134\,890) "
               r"graphs, scored on the nodes of their ground-truth communities (top 5000 for Amazon and DBLP, "
               r"all groups of at least three nodes for YouTube), to which every detected community is "
               r"restricted. Best-match F1 against the recorded communities and $\Omega$ against their distinct "
               r"node sets, at each method's best parameter by that metric (best) and at the parameter chosen "
               r"without labels by maximal overlapping modularity on the full graph (free). HC$_f$: the seven "
               r"commonality functions together, the free choice taken over all of them. Connected components "
               r"are a trivial reference. On YouTube HyperCommon uses degree-ratio pruning (exact for "
               r"Jaccard), which excludes link communities and the other functions; clique percolation and SLPA "
               r"are compared on the induced subgraphs (Table~\ref{tab:realinduced}). Bold: best in column.")
        label, name = "tab:real", "table_real.tex"
    else:
        cap = (r"As Table~\ref{tab:real}, detecting on the subgraphs induced by the labelled nodes (16\,716 "
               r"nodes for Amazon, 93\,432 for DBLP, 47\,135 for YouTube). The induced Amazon subgraph falls into "
               r"1\,106 components that already score F1 $0.958$, and on YouTube the components are the most "
               r"accurate cover by F1, so this setting does not discriminate between methods on those two. "
               r"Clique percolation on YouTube did not complete (memory above 100\,GB). Bold: best in column "
               r"among the methods; the components row is a reference and is not marked.")
        label, name = "tab:realinduced", "table_real_induced.tex"
    out = [r"\begin{table*}[t]", r"\caption{" + cap + "}", r"\label{" + label + "}",
           r"\scriptsize\setlength{\tabcolsep}{2.5pt}", r"\begin{ruledtabular}",
           r"\begin{tabular}{l" + "c" * 12 + "}",
           " & " + " & ".join(r"\multicolumn{4}{c}{" + n + "}" for _, n in REAL_NETS) + r" \\",
           " & " + " & ".join(["F1 best & F1 free & $\\Omega$ best & $\\Omega$ free"] * 3) + r" \\", r"\hline"]
    for m, lab, vals in rows:
        cc = ["--" if np.isnan(v) else (r"\textbf{%.3f}" % v if m != "components" and abs(v - best[j]) < 5e-4
                                          else "%.3f" % v) for j, v in enumerate(vals)]
        out.append(lab + " & " + " & ".join(cc) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    write(name, out)


MORE_NETS = [("football", "College football"), ("polbooks", "Political books"), ("polblogs", "Political blogs"),
             ("primary_school", "Primary school")] + [
    (f"facebook_{e}", f"Facebook ego {e}") for e in (0, 107, 348, 414, 686, 698, 1684, 1912, 3437, 3980)]
MORE_OTHERS = ["linkcomm", "cpm", "leiden", "leiden_cpm", "infomap", "lpa", "walktrap", "slpa", "demon", "angel",
               "scan", "scan_j"]


def more_picks(metric):
    """Per network and method: the metric at the best parameter (mean over seeds of the per-seed
    best) and at the label-free parameter (per seed, largest eq); HC_f over the seven functions."""
    raw = pd.read_csv(os.path.join(RES, "real_more_v3.csv"), dtype={"param": str})
    d = raw[raw.error.isna()].copy()
    fam = d[d.method.str.startswith("hcf_")].copy()
    fam["method"] = "hc_family"
    d = pd.concat([d, fam], ignore_index=True)
    key = ["network", "method", "seed"]
    best = d.loc[d.groupby(key)[metric].idxmax()].groupby(["network", "method"])[metric].mean().unstack()
    free = d.loc[d.groupby(key)["eq"].idxmax()].groupby(["network", "method"])[metric].mean().unstack()
    return best, free, raw[raw.error.notna()][["network", "method"]].drop_duplicates()


def more_table():
    out_all = []
    for metric, label, mname in (("f1", "tab:realmore", "best-match F1"), ("omega", "tab:realmoreomega", r"$\Omega$")):
        best, free, failed = more_picks(metric)
        fails = ", ".join(f"{NAME.get(m, m)} on {dict(MORE_NETS)[n]}" for n, m in failed.itertuples(index=False))
        cap = (r"Detection on fourteen further networks with recorded groups (Netzschleuder "
               r"repository~\cite{peixoto2020netzschleuder}; sizes in Table~\ref{tab:datamore}), scored by "
               + mname + r" on the labelled nodes, to which detected communities are restricted; each method at "
               r"its best parameter (best) and at the parameter of largest overlapping modularity (free); "
               r"stochastic methods: mean over five seeds. HC: HyperCommon with Jaccard; HC$_f$: the seven "
               r"commonality functions together, each on 100 quantile levels of its own triple scores (in "
               r"brackets the function chosen; J when Jaccard ties), the free choice taken over all of them; with "
               r"labels it chooses among 700 (function, threshold) pairs per network, so its best-parameter "
               r"column is an upper bound for the family, not one method's parameter grid (with 30 quantile "
               r"levels per function its four wins are unchanged). "
               r"``Best other'': the most accurate of the twelve other methods under that protocol"
               + (f" (runs that did not complete: {fails})" if fails and metric == "f1" else "") +
               r". Bold: best of HC, HC$_f$ and the other methods under each protocol. Functions: J, Jaccard; "
               r"cos, cosine; $I/\min$ and $I/\max$, overlap normalised by the smaller or larger neighbourhood; $I$, overlap "
               r"count; RA, resource allocation. College football: the conference labels of "
               r"Ref.~\cite{girvan2002community} are those of the 2001 season while the games are from the 2000 "
               r"season; corrected labels are given in Ref.~\cite{evans2010clique}. Every method on every network "
               r"is in Tables~\ref{tab:S-realnet-f1-best}--\ref{tab:S-realnet-coverage-free}.")
        if metric == "omega":
            cap = r"As Table~\ref{tab:realmore}, scored by $\Omega$ on the labelled nodes."
        out = [r"\begin{table*}[t]", r"\caption{" + cap + "}", r"\label{" + label + "}", r"\footnotesize",
               r"\begin{ruledtabular}", r"\begin{tabular}{l|ccl|ccl}",
               r" & \multicolumn{3}{c|}{best parameter} & \multicolumn{3}{c}{free (label-free)} \\",
               r"Network & HC & HC$_f$ & best other & HC & HC$_f$ & best other \\", r"\hline"]
        raw = pd.read_csv(os.path.join(RES, "real_more_v3.csv"), dtype={"param": str})
        raw = raw[raw.error.isna() & raw.method.str.startswith("hcf_")]
        rows = []
        for key, lab in MORE_NETS:
            line, row = lab, dict(network=key)
            for prot, T in (("best", best), ("free", free)):
                r = T.loc[key]
                oth = r[[m for m in MORE_OTHERS if m in r.index]].dropna()
                ob, on = oth.max(), oth.idxmax()
                hc, hf = r["hc"], r["hc_family"]
                g = raw[raw.network == key]
                if prot == "best":
                    pf = g.loc[g.groupby("seed")[metric].idxmax()].method.mode().iloc[0][4:]
                else:
                    pf = g.loc[g.groupby("seed")["eq"].idxmax()].method.mode().iloc[0][4:]
                fname = "J" if abs(hf - r.get("hcf_jaccard", -1)) < 5e-4 else F_SHORT[pf]
                top = max(hc, hf, ob)
                fmt = lambda v: ("$-$%.3f" % -v) if v < 0 else "%.3f" % v
                b = lambda v: (r"\textbf{" + fmt(v) + "}") if abs(v - top) < 5e-4 else fmt(v)
                line += f" & {b(hc)} & {b(hf)} ({fname}) & {b(ob)} ({NAME.get(on, on)})"
                row.update({f"hc_{prot}": hc, f"hcf_{prot}": hf, f"f_{prot}": pf, f"other_{prot}": ob,
                            f"other_method_{prot}": on})
            out.append(line + r" \\")
            rows.append(row)
        out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
        pd.DataFrame(rows).to_csv(os.path.join(DATA, f"real_more_summary2_{metric}.csv"), index=False)
        out_all += out + [""]
    write("table_real_more.tex", out_all)


F_SHORT = {"jaccard": "J", "cosine": "cos", "max_norm": r"$I/\max$", "min_norm": r"$I/\min$", "i2_over_u": r"$I^2/U$",
           "count": "$I$", "resource_alloc": "RA"}

FAM_FUNCS = [("jaccard", r"Jaccard"), ("cosine", "cosine"), ("max_norm", r"$I/\max$"),
             ("min_norm", r"$I/\min$"), ("i2_over_u", r"$I^2/U$"), ("count", r"$I$"),
             ("resource_alloc", "RA")]


def family():
    """HyperCommon with each commonality function (results2/family2.csv): mean Omega over the
    noise range under the oracle and label-free protocols, paired difference from Jaccard."""
    sys.path.insert(0, BENCH)
    from analysis2 import select, areas, boot
    parts = [pd.read_csv(os.path.join(R2, n), dtype={"param": str}) for n in ("family2.csv", "family2_top.csv")
             if os.path.exists(os.path.join(R2, n))]
    df = pd.concat(parts, ignore_index=True)
    df = df[df.error.isna()].copy()
    df["param"] = df["param"].fillna("")
    sel, _ = select(df)
    anyf, _ = select(df.assign(method="hcf_any"))   # function and level chosen together per graph
    label = {f: lab for _, f, lab in ROWS}
    fams = [f for _, f, _ in ROWS if f in set(df.family)]
    rows, out_all = [], []
    for prot, tag, word in (("oracle", "familyoracle", "best threshold per graph"),
                            ("label_free", "familyfree", "label-free threshold (largest overlapping modularity)")):
        g = sel[sel.protocol == prot]
        ga = anyf[anyf.protocol == prot]
        out = [r"\begin{table*}[t]",
               r"\caption{HyperCommon with each commonality function, " + word + ": mean $\\Omega$ over the "
               r"noise range ($p\in[0,0.8]$ for rings, $\mu\in[0.1,0.6]$ for LFR), ten test graphs per noise "
               r"level, the threshold a quantile of each function's triple scores. Superscripts mark a paired "
               r"95\% interval of the difference from Jaccard that excludes zero ($^+$ better, $^-$ worse). "
               r"Bold: best single function in the row (last column excluded). $I/\max$, $I/\min$: overlap "
               r"normalised by the larger or smaller neighbourhood; $I$: overlap count; RA: resource allocation. "
               r"Last column: function and level chosen together on each graph by the same rule. Quantile levels "
               r"$0,0.04,\dots,0.96$ and $0.97,0.98,0.99,0.995$; rings sampled every $0.1$ in $p$. Jaccard values "
               r"can differ from the HC columns of Tables~\ref{tab:labelfree} and~\ref{tab:oracle} (by up to "
               r"$0.04$) because thresholds are quantile levels and rings are sampled every $0.1$ in $p$.}",
               r"\label{tab:" + tag + "}", r"\scriptsize", r"\begin{ruledtabular}",
               r"\begin{tabular}{l" + "c" * len(FAM_FUNCS) + "|c}",
               "Family & " + " & ".join(lab for _, lab in FAM_FUNCS) + r" & per graph \\", r"\hline"]
        wins = {f: [0, 0] for f, _ in FAM_FUNCS}
        for fam in fams:
            a = {f: areas(g[(g.family == fam) & (g.method == "hcf_" + f)], "omega") for f, _ in FAM_FUNCS}
            ref = a["jaccard"]
            vals, sig = {}, {}
            for f, _ in FAM_FUNCS:
                vals[f] = a[f].mean()
                if f == "jaccard":
                    continue
                c = a[f].index.intersection(ref.index)
                m, lo, hi = boot((a[f][c] - ref[c]).values)
                sig[f] = "+" if lo > 0 else "-" if hi < 0 else ""
                wins[f][0] += sig[f] == "+"
                wins[f][1] += sig[f] == "-"
                rows.append(dict(protocol=prot, family=fam, function=f, mean=vals[f], diff=m, lo=lo, hi=hi))
            rows.append(dict(protocol=prot, family=fam, function="jaccard", mean=vals["jaccard"]))
            best = max(vals.values())
            cells = []
            for f, _ in FAM_FUNCS:
                s = f"{vals[f]:.2f}" + (f"$^{sig[f]}$" if sig.get(f) else "")
                cells.append(r"\textbf{" + s + "}" if round(vals[f], 2) == round(best, 2) else s)
            va = areas(ga[ga.family == fam], "omega").mean()
            rows.append(dict(protocol=prot, family=fam, function="per_graph", mean=va))
            out.append(label[fam] + " & " + " & ".join(cells) + f" & {va:.2f}" + r" \\")
        out.append(r"\hline")
        out.append("better / worse than Jaccard & -- & " +
                   " & ".join(f"{wins[f][0]} / {wins[f][1]}" for f, _ in FAM_FUNCS[1:]) + r" & \\")
        out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}", ""]
        out_all += out
    write("table_family.tex", out_all)
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "family2_vs_jaccard.csv"), index=False)


if __name__ == "__main__":
    which = sys.argv[1:] or ["benchmarks", "robust", "real", "more"]
    if "benchmarks" in which:
        benchmarks()
    if "robust" in which:
        robust()
    if "real" in which:
        real_table("full")
        real_table("induced")
    if "more" in which:
        more_table()
    if "family" in which:
        family()
    if "geometry" in which:
        geometry()
    print("done:", which)
