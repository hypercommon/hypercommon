"""
Generate the results tables of the paper straight from the benchmark summaries,
so no number in a table is typed by hand.

  sections/table_benchmarks.tex   mean Omega over the noise range, oracle protocol,
                                  for every benchmark family x selected methods
  figures/data/table_benchmarks.csv   the same numbers, with 95% CIs
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "bench", "results")
SEC = os.path.join(HERE, "..", "sections")
DATA = os.path.join(HERE, "..", "figures", "data")

FAMILIES = [  # (experiment, family, label)
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
METHODS = [("hc", "HC"), ("hc_tri", "tri."), ("hc_2test", "2-test"), ("cpm", "CPM"),
           ("leiden", "Leiden"), ("infomap", "Infomap"), ("slpa", "SLPA"), ("angel", "ANGEL"),
           ("demon", "DEMON"), ("scan_j", "SCAN"), ("linkcomm", "Link")]


def load(protocol, metric="omega"):
    if metric == "onmi":
        rows = []
        for exp in {f[0] for f in FAMILIES}:
            p = os.path.join(RES, f"{exp}_summary.json")
            if not os.path.exists(p):
                continue
            for x in json.load(open(p))["summary"]:
                if x["protocol"] != protocol:
                    continue
                ps = np.array(x["p"])
                rows.append(dict(experiment=exp, family=x["family"], method=x["method"],
                                 mean_omega=float(np.trapezoid(x["onmi"], ps) / (ps[-1] - ps[0])), lo=np.nan, hi=np.nan))
        return pd.DataFrame(rows)
    rows = []
    for exp in {f[0] for f in FAMILIES}:
        p = os.path.join(RES, f"{exp}_summary.json")
        if not os.path.exists(p):
            continue
        for x in json.load(open(p))["summary"]:
            if x["protocol"] != protocol:
                continue
            span = x["p"][-1] - x["p"][0]
            rows.append(dict(experiment=exp, family=x["family"], method=x["method"],
                             mean_omega=x["area"] / span, lo=x["area_lo"] / span, hi=x["area_hi"] / span))
    return pd.DataFrame(rows)


def table(protocol, label, caption, metric="omega"):
    df = load(protocol, metric)
    df.to_csv(os.path.join(DATA, f"table_benchmarks_{metric}_{protocol}.csv"), index=False)
    cols = "l" + "c" * len(METHODS)
    out = [r"\begin{table*}[t]", r"\caption{" + caption + "}", r"\label{" + label + "}",
           r"\footnotesize\setlength{\tabcolsep}{2.5pt}",
           r"\begin{ruledtabular}", r"\begin{tabular}{" + cols + "}",
           "Family & " + " & ".join(m[1] for m in METHODS) + r" \\", r"\hline"]
    last_exp = None
    for exp, fam, lab in FAMILIES:
        sub = df[(df.experiment == exp) & (df.family == fam)].set_index("method")
        if sub.empty:
            continue
        if last_exp is not None and exp != last_exp and not (last_exp.startswith("lfr") and exp.startswith("lfr")):
            out.append(r"\hline")
        last_exp = exp
        vals = {m: sub.mean_omega.get(m, np.nan) for m, _ in METHODS}
        best = np.nanmax(list(vals.values()))
        cells = []
        for m, _ in METHODS:
            v = vals[m]
            if np.isnan(v):
                cells.append("--")
            else:
                s = f"{v:.2f}"
                cells.append(r"\textbf{" + s + "}" if abs(v - best) < 0.005 else s)
        out.append(lab + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    return "\n".join(out)


if __name__ == "__main__":
    d = load("oracle")
    hw = (d.hi - d.lo) / 2
    ci_note = (rf"95\% bootstrap half-widths have median {hw.median():.3f} and maximum {hw.max():.3f} "
               r"(all intervals in the supplementary data).")
    tex = table("oracle", "tab:oracle",
                r"Mean $\Omega$ over the noise range ($p\in[0,0.8]$ for rings, $\mu\in[0.1,0.6]$ for "
                r"LFR, $\mu\in[0.2,0.6]$ for the LFR slices), each method at its best parameter on each "
                r"test graph (oracle protocol), ten test graphs per noise level. HC: HyperCommon; tri.: "
                r"triangles only (Jaccard-filtered 3-clique percolation); 2-test: third pair not tested; "
                r"SCAN with Jaccard similarity. Bold: best in row, within 0.005. " + ci_note)
    tex += "\n\n" + table("frozen", "tab:frozen",
                          r"As Table~\ref{tab:oracle}, with each method's parameter frozen on calibration "
                          r"graphs of the same family, all noise levels pooled (frozen protocol).")
    open(os.path.join(SEC, "table_benchmarks.tex"), "w", encoding="utf-8").write(tex + "\n")
    tex_onmi = table("oracle", "tab:onmi",
                     r"As Table~\ref{tab:oracle}, scored with ONMI (McDaid \emph{et al.}, max normalisation) "
                     r"on the same covers, i.e.\ at each method's $\Omega$-optimal parameter.", metric="onmi")
    open(os.path.join(SEC, "table_onmi.tex"), "w", encoding="utf-8").write(tex_onmi + "\n")
    print(tex[:2500])


FAMILY_COLS = [("rings_k8_o05", r"rings $k{=}8$"), ("rings_k16_o05", r"rings $k{=}16$"),
               ("two_worlds", "two worlds"), ("gradient10", "gradient"),
               ("lfr_on10", r"LFR $o_n{=}0.1$"), ("lfr_on30", r"LFR $o_n{=}0.3$"),
               ("lfr_k10", r"LFR $\langle k\rangle{=}10$")]
FAMILY_ROWS = [("jaccard", r"Jaccard $I/U$"), ("cosine", r"cosine $I/\sqrt{|N[u]||N[v]|}$"),
               ("max_norm", r"$I/\max(|N[u]|,|N[v]|)$"), ("min_norm", r"$I/\min(|N[u]|,|N[v]|)$"),
               ("i2_over_u", r"$I^2/U$ (and $I/U^{1/2}$)"), ("count", r"overlap $I$"),
               ("resource_alloc", "resource allocation")]


def family_table():
    S = json.load(open(os.path.join(RES, "commonality_family_summary.json")))
    by = {(x["family"], x["method"][4:]): x for x in S["summary"] if x["protocol"] == "oracle"}
    pair = pd.read_csv(os.path.join(DATA, "family_vs_jaccard_paired.csv"))
    sig = {(r.family, r.method): ("+" if r.lo > 0 else "-" if r.hi < 0 else "") for r in pair.itertuples()}
    mean = {k: v["area"] / (v["p"][-1] - v["p"][0]) for k, v in by.items()}
    out = [r"\begin{table*}[t]",
           r"\caption{HyperCommon with each commonality function: mean $\Omega$ over the noise range "
           r"($p\in[0,0.6]$ for rings, $\mu\in[0.1,0.6]$ for LFR), best threshold per graph, the threshold "
           r"taken as a quantile of each function's triple scores. Superscripts mark a paired 95\% interval "
           r"of the difference from Jaccard that excludes zero ($^+$ better, $^-$ worse). Bold: best in "
           r"column. $I/U^{1/2}=(I^2/U)^{1/2}$ ranks pairs identically to $I^2/U$.}",
           r"\label{tab:family}", r"\footnotesize", r"\begin{ruledtabular}",
           r"\begin{tabular}{l" + "c" * len(FAMILY_COLS) + "}",
           "Function & " + " & ".join(c[1] for c in FAMILY_COLS) + r" \\", r"\hline"]
    best = {f: max(mean[(f, m)] for m, _ in FAMILY_ROWS) for f, _ in FAMILY_COLS}
    for m, lab in FAMILY_ROWS:
        cells = []
        for f, _ in FAMILY_COLS:
            v = mean[(f, m)]
            s = f"{v:.3f}" + (f"$^{sig[(f, m)]}$" if sig.get((f, m)) else "")
            cells.append(r"\textbf{" + s + "}" if abs(v - best[f]) < 5e-4 else s)
        out.append(lab + " & " + " & ".join(cells) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    open(os.path.join(SEC, "table_family.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    pd.DataFrame([dict(family=f, function=m, mean_omega=mean[(f, m)]) for (f, m) in mean]).to_csv(
        os.path.join(DATA, "family_mean_omega.csv"), index=False)
    return "\n".join(out)


MORE_NETS = [("football", "College football"), ("polbooks", "Political books"), ("polblogs", "Political blogs"),
             ("primary_school", "Primary school")] + [
    (f"facebook_{e}", f"Facebook ego {e}") for e in (0, 107, 348, 414, 686, 698, 1684, 1912, 3437, 3980)]
F_SHORT = {"jaccard": "J", "cosine": "cos", "max_norm": "max", "min_norm": "min", "i2_over_u": r"$I^2/U$",
           "i_over_sqrtu": r"$I^2/U$", "count": "$I$", "resource_alloc": "RA"}
OTHERS = ["cpm", "leiden", "infomap", "lpa", "walktrap", "slpa", "demon", "angel", "scan", "scan_j", "linkcomm"]
OTHER_NAME = {"cpm": "CPM", "leiden": "Leiden", "infomap": "Infomap", "lpa": "LPA", "walktrap": "Walktrap",
              "slpa": "SLPA", "demon": "DEMON", "angel": "ANGEL", "scan": "SCAN", "scan_j": "SCAN-J",
              "linkcomm": "Link"}


def real_more_table():
    import networkx as nx
    import sys as _sys
    _sys.path.insert(0, HERE)
    from real_more import load
    raw = pd.read_csv(os.path.join(RES, "real_more_v2.csv"), dtype={"error": str})
    failed = raw[raw.error.notna()][["network", "method", "error"]].drop_duplicates(["network", "method"])
    d = raw[raw.error.isna()]
    rules = {}
    for col in ("f1", "f1_ge3"):
        best = d.loc[d.groupby(["network", "method", "seed"])[col].idxmax()].groupby(["network", "method"])[col].mean()
        rules[col] = best.unstack("method")
    hcf = [c for c in rules["f1"].columns if c.startswith("hcf_")]
    # runs that raised MemoryError, and runs stopped by hand once they had exceeded physical memory
    stopped = failed.error.astype(str).str.startswith("MemoryError: stopped by operator")
    name = lambda f: ", ".join(f"{OTHER_NAME.get(m, m)} on {dict(MORE_NETS)[n]}"
                               for n, m in f[["network", "method"]].itertuples(index=False))
    fails = "; ".join(s for s in (
        (f"out of memory: {name(failed[~stopped])}" if (~stopped).any() else ""),
        (f"stopped once together they had exhausted the machine's 31\\,GB of memory:{name(failed[stopped])}" if stopped.any() else ""))
        if s)
    out = [r"\begin{table*}[t]",
           r"\caption{Detection on fourteen further networks with recorded groups (Netzschleuder "
           r"repository~\cite{peixoto2020netzschleuder}; sizes "
           r"in Table~\ref{tab:datamore}). Best-match F1 on the labelled nodes, to which detected communities are "
           r"restricted, keeping every nonempty restricted community (left) or only those of at least three nodes "
           r"(right); each method at its best parameter under each rule (stochastic methods: mean over five seeds "
           r"of the per-seed best). Every HyperCommon variant searches the same 100 quantile levels of its own "
           r"triple scores. comp.: the connected components returned as communities. HC: HyperCommon with "
           r"Jaccard; HC$_f$: with the best member of the commonality family on that network (in brackets; J "
           r"when Jaccard ties) -- a selection made with the labels, reported as what the family can reach, not "
           r"as a method. ``Best other'': the most accurate of the eleven other methods" +
           (f"; runs that did not complete -- {fails}" if fails else "") +
           r". $C$: transitivity. Bold: best of HC and the other methods under each rule. On political blogs "
           r"every HyperCommon variant does best with a single community, the giant component.}",
           r"\label{tab:realmore}", r"\footnotesize", r"\begin{ruledtabular}",
           r"\begin{tabular}{lcr|cccl|cccl}",
           r" & & & \multicolumn{4}{c|}{all restricted communities} & "
           r"\multicolumn{4}{c}{restricted communities of $\ge3$ nodes} \\",
           r"Network & $C$ & groups & comp. & HC & HC$_f$ & best other & comp. & HC & HC$_f$ & best other \\",
           r"\hline"]
    rows, sizes = [], []
    for key, lab in MORE_NETS:
        if key not in rules["f1"].index:
            continue
        G, truth = load(key)
        line = f"{lab} & {nx.transitivity(G):.2f} & {len(truth)}"
        row = dict(network=key)
        for col, F in rules.items():
            r = F.loc[key]
            hc = r["hc"]
            fbest = r[hcf].max()
            # on a tie with Jaccard, name Jaccard
            fname = "J" if abs(r["hcf_jaccard"] - fbest) < 5e-4 else F_SHORT[r[hcf].idxmax()[4:]]
            oth = r[[m for m in OTHERS if m in r.index]].dropna()
            ob, on = oth.max(), OTHER_NAME[oth.idxmax()]
            top = max(hc, ob)
            cell = lambda v: (r"\textbf{" + f"{v:.3f}" + "}") if abs(v - top) < 5e-4 else f"{v:.3f}"
            line += f" & {r['components']:.3f} & {cell(hc)} & {fbest:.3f} ({fname}) & {cell(ob)} ({on})"
            sfx = "" if col == "f1" else "_ge3"
            row.update({f"hc{sfx}": hc, f"hc_best_f{sfx}": fbest, f"best_f{sfx}": r[hcf].idxmax()[4:],
                        f"best_other{sfx}": ob, f"best_other_method{sfx}": oth.idxmax()})
        out.append(line + r" \\")
        rows.append(row)
        lab_nodes = set().union(*truth)
        sizes.append(f"{lab} & {G.number_of_nodes()} & {G.number_of_edges()} & {len(truth)} & "
                     f"{len(lab_nodes & set(G)) / G.number_of_nodes():.2f} \\\\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    open(os.path.join(SEC, "table_real_more.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "real_more_summary.csv"), index=False)
    data = [r"\begin{table}[h]",
            r"\caption{The fourteen networks of Table~\ref{tab:realmore}: nodes, edges, recorded groups, and "
            r"the fraction of nodes that carry at least one group; groups of fewer than three nodes are "
            r"discarded. Political blogs includes 266 isolated nodes; the primary school network aggregates "
            r"both days of contacts.}",
            r"\label{tab:datamore}", r"\begin{ruledtabular}", r"\begin{tabular}{lrrrr}",
            r"Network & Nodes & Edges & Groups & $f_{\rm lab}$ \\", r"\hline"] + sizes + [
            r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    open(os.path.join(SEC, "table_data_more.tex"), "w", encoding="utf-8").write("\n".join(data) + "\n")
    return "\n".join(out)


def overlap_table():
    d = pd.read_csv(os.path.join(DATA, "hetero_by_overlap.csv"))
    names = {"two_worlds": "two worlds", "three_big_many_small": "three big, many small",
             "gradient10": "gradient of ten rings", "sizes_only_k8": r"equal $k$, unequal $s$ ($k{=}8$)"}
    mark = {"ahead": "", "tie": r"$^{\approx}$", "behind": r"$^-$"}
    out = [r"\begin{table*}[t]",
           r"\caption{Unequal rings at three overlap levels: mean $\Omega$ over $p\in[0,0.8]$, best parameter per "
           r"graph, for HyperCommon (HC), the most accurate other method, and Leiden. Marks on HC: $^{\approx}$ "
           r"paired 95\% interval of the difference from the best other method includes zero, $^-$ HC behind; "
           r"unmarked, HC ahead.}", r"\label{tab:overlaplevels}", r"\footnotesize\setlength{\tabcolsep}{3pt}",
           r"\begin{ruledtabular}", r"\begin{tabular}{lcccl}",
           r"Family & $o_n$ & HC & Leiden & best other \\", r"\hline"]
    for fam in ("two_worlds", "three_big_many_small", "gradient10", "sizes_only_k8"):
        for r in d[d.family == fam].sort_values("overlap").itertuples():
            out.append(f"{names[fam] if r.overlap == 0 else ''} & {r.overlap / 100:.2f} & {r.hc:.3f}{mark[r.vs_best]} & "
                       f"{r.leiden:.3f} & {r.best_other_val:.3f} ({OTHER_NAME.get(r.best_other, r.best_other)}) \\\\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    open(os.path.join(SEC, "table_overlap.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    return "\n".join(out)


REAL_CAPTION = {
    "full": (r"Detection on real networks, run on the full SNAP Amazon (334\,863 nodes) and DBLP "
             r"(317\,080 nodes) graphs and scored on the nodes of their top 5000 ground-truth communities, "
             r"to which every detected community is restricted; every nonempty restricted community is kept. "
             r"Best-match F1, ONMI and coverage at each method's best parameter by F1 (grids widened until the "
             r"best value was interior); F1$_{\ge3}$: best F1 when restricted communities of fewer than three "
             r"nodes are also dropped, a sensitivity check; ``free'' is F1 without any "
             r"label: HyperCommon and the triangles-only variant at the threshold maximising overlapping "
             r"modularity, Leiden at $\gamma=1$, parameter-free methods as they are. Connected components are "
             r"a trivial reference. Clique percolation and SLPA, too slow on the full graphs in their "
             r"available implementations, are compared on the induced subnetworks only "
             r"(Table~\ref{tab:realinduced}). Bold: best in column.", "tab:real"),
    "induced": (r"As Table~\ref{tab:real}, but detecting on the subnetwork induced by the ground-truth "
                r"communities (16\,716 nodes for Amazon, 93\,432 for DBLP), as in "
                r"Ref.~\cite{yang2015defining}. The induced Amazon subnetwork falls into 1\,106 components that "
                r"already score F1 $0.958$, so this setting does not discriminate between methods on Amazon.",
                "tab:realinduced"),
}


def real_table(setting="full"):
    """Per network: F1 and ONMI at the best parameter by F1, coverage of that cover, F1 at
    the best parameter under the >=3 rule (sensitivity), and F1 without labels."""
    d = pd.read_csv(os.path.join(RES, "real_detection_v3.csv"))
    d = d[d.setting == setting]
    names = [("components", "components"), ("hc", "HyperCommon"), ("hc_tri", "triangles only"),
             ("cpm", "CPM"), ("leiden", "Leiden"), ("infomap", "Infomap"), ("slpa", "SLPA"),
             ("lpa", "label prop.")]
    K = 5
    rows = []
    for m, lab in names:
        cells = []
        for net in ("amazon", "dblp"):
            g = d[(d.network == net) & (d.method == m)]
            if g.empty:
                cells += ["--"] * K
                continue
            b = g.loc[g.f1.idxmax()]
            lf = g[g.label_free.notna() & (g.label_free != "")]
            lf_f1 = f"{lf.f1.iloc[0]:.3f}" if len(lf) else "--"
            cells += [f"{b.f1:.3f}", f"{b.onmi:.3f}" if pd.notna(b.onmi) else "--", f"{b.coverage:.2f}",
                      f"{g.f1_ge3.max():.3f}", lf_f1]
        if any(c != "--" for c in cells):    # methods not run in this setting are named in the caption
            rows.append((lab, cells))
    best = {}
    for j in range(2 * K):
        if j % K == 2:                      # no bold for coverage
            best[j] = None
            continue
        vals = [float(c[j]) for lab, c in rows if c[j] != "--" and lab != "components"]
        best[j] = max(vals) if vals else None
    cap, lab_ = REAL_CAPTION[setting]
    out = [r"\begin{table*}[t]", r"\caption{" + cap + "}",
           r"\label{" + lab_ + "}", r"\footnotesize",
           r"\begin{ruledtabular}", r"\begin{tabular}{l" + "c" * (2 * K) + "}",
           r" & \multicolumn{5}{c}{Amazon} & \multicolumn{5}{c}{DBLP} \\",
           r"Method & F1 & ONMI & cov. & F1$_{\ge3}$ & free & F1 & ONMI & cov. & F1$_{\ge3}$ & free \\", r"\hline"]
    for lab, cells in rows:
        cc = [(r"\textbf{" + c + "}" if c != "--" and best[j] is not None and abs(float(c) - best[j]) < 5e-4 else c)
              for j, c in enumerate(cells)]
        out.append(lab + " & " + " & ".join(cc) + r" \\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    name = "table_real.tex" if setting == "full" else "table_real_induced.tex"
    open(os.path.join(SEC, name), "w", encoding="utf-8").write("\n".join(out) + "\n")
    return "\n".join(out)


if __name__ == "__main__":
    print(real_table("full"))
    print(real_table("induced"))
    print(family_table())
    print(real_more_table())
    print(overlap_table())


def double_bond_table():
    """Appendix: ring pairs sharing a close pair of merged nodes (figures/data/double_bonds.csv,
    analysis/double_bonds.py)."""
    d = pd.read_csv(os.path.join(DATA, "double_bonds.csv"))
    names = {"two_worlds": "two worlds", "three_big_many_small": "three big, many small",
             "gradient10": "gradient", "dense_even": "equal $k$, unequal $s$ ($k{=}16$)",
             "sizes_only_k8": "equal $k$, unequal $s$ ($k{=}8$)", "k_only_s100": "equal $s$, unequal $k$"}
    out = [r"\begin{table}[t]",
           r"\caption{Shared pairs in the ring benchmarks without noise (mean over ten graphs): merges, ring pairs "
           r"sharing at least one node, ring pairs sharing two merged nodes within ring distance $k$ on both "
           r"rings (Proposition~\ref{prop:sharedpairs}), the share of rings in such a pair, and the share of rings "
           r"that exact HyperCommon fuses with another ring just above the lower edge $1/(2k_{\min}+1)$ of the "
           r"recovery window (a detected community holding at least three unshared nodes of each of two rings; "
           r"analysis/double\_bonds.py).}",
           r"\label{tab:doublebonds}", r"\scriptsize\setlength{\tabcolsep}{2pt}",
           r"\begin{ruledtabular}", r"\begin{tabular}{lcrrrrr}",
           r"Family & $o_n$ & merges & sharing & pair & rings & fused \\", r"\hline"]
    for r in d.itertuples():
        fam = r.family
        o = 0.10 if fam.endswith("_o10") or fam.endswith("o10") else 0.05
        base = fam.replace("_o10", "").replace("_o05", "")
        if base.startswith("k8_s50"):
            lab = r"uniform $(50,8)$"
        elif base.startswith("k16_s100"):
            lab = r"uniform $(100,16)$"
        else:
            lab = names.get(base, base)
        out.append(f"{lab} & {o:.2f} & {r.merges:.0f} & {r.ring_pairs_sharing:.1f} & "
                   f"{r.double_bonded_pairs:.1f} & {r.share_rings_double_bonded:.2f} & {r.share_rings_fused:.2f} \\\\")
    out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    open(os.path.join(SEC, "table_doublebonds.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    return "\n".join(out)
