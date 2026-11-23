"""
Control asked for in review: does the triple construction add anything beyond the commonality
itself? Leiden-RB modularity on the graph weighted by the closed-neighbourhood Jaccard of each
edge (leiden_j, 30 resolutions as Leiden; bench/methods2.py, rows in
bench/results2/<exp>_jctrl.csv from queue_review.sh), against HyperCommon and unweighted Leiden,
on the 28 ring settings and the eight LFR families of the scorecard. Mean Omega over the noise
range and the ten test graphs per noise level (as the scorecard), oracle (best parameter per
graph) and label-free (largest overlapping modularity); HC minus leiden_j with a paired 95%
bootstrap interval over the ten test seeds, each averaged over the noise range (as the scorecard).

    python jaccard_control.py
Writes figures/data/jaccard_control.csv and sections/table_jaccard_control.tex.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.join(HERE, "..", "bench", "results2")
DATA = os.path.join(HERE, "..", "figures", "data")
TEX = os.path.join(HERE, "..", "sections", "table_jaccard_control.tex")
rng = np.random.default_rng(11)
COLS = ["family", "p", "seed", "split", "method", "omega", "eq"]
RING_EXPS = ["rings_uniform", "rings_hetero", "rings_hetero_overlap", "rings_pref"]
LFR_EXPS = ["lfr", "lfr_slices", "lfr_clustered"]
LFR_FAMS = {("lfr", "on00"): r"LFR, $o_n{=}0$", ("lfr", "on10"): r"LFR, $o_n{=}0.1$",
            ("lfr", "on30"): r"LFR, $o_n{=}0.3$", ("lfr_slices", "k10"): r"LFR, $\langle k\rangle{=}10$",
            ("lfr_slices", "om4"): r"LFR, $o_m{=}4$", ("lfr_slices", "big"): r"LFR, sizes 50--250",
            ("lfr_clustered", "on10_C03"): r"LFR, $C{=}0.3$", ("lfr_clustered", "on10_C05"): r"LFR, $C{=}0.5$"}
BASE = {"k8_s50": r"uniform, $k{=}8$", "k16_s100": r"uniform, $k{=}16$",
        "sizes_only_k8": r"equal $k$, unequal $s$, $k{=}8$", "dense_even": r"equal $k$, unequal $s$, $k{=}16$",
        "k_only_s100": r"equal $s$, unequal $k$", "two_worlds": "two worlds",
        "three_big_many_small": "three big, many small", "gradient10": "gradient"}
ORDER = list(BASE)


def ring_name(exp, fam):
    if exp == "rings_pref":
        b = BASE["k8_s50"] if fam.startswith("k8") else BASE["two_worlds"]
        return b + (r", $\alpha{=}2$" if fam.endswith("2") else r", $\alpha{=}1$"), "0.05", (1, 0 if fam.startswith("k8") else 1, fam)
    key = next(k for k in ORDER if fam.startswith(k))
    rest = fam[len(key):]
    on = {"_o00": "0", "_o05": "0.05", "_o10": "0.10"}.get(rest[:4], "0.05")
    return BASE[key], on, (0, ORDER.index(key), on)


def load(exp):
    main = pd.read_csv(os.path.join(R2, f"{exp}.csv"), usecols=COLS, low_memory=False)
    main = main[main.method.isin(["hc", "leiden"])]
    ctrl = pd.read_csv(os.path.join(R2, f"{exp}_jctrl.csv"), usecols=COLS, low_memory=False)
    d = pd.concat([main, ctrl], ignore_index=True)
    return d[(d.split == "test") & d.omega.notna()]


def selected(d):
    out = {}
    for prot, col in (("oracle", "omega"), ("label_free", "eq")):
        x = d[d[col].notna()]
        s = x.loc[x.groupby(["family", "method", "p", "seed"])[col].idxmax()]
        out[prot] = s.pivot_table(index=["family", "p", "seed"], columns="method", values="omega")
    return out


def main():
    rows = []
    for exp in RING_EXPS + LFR_EXPS:
        sel = selected(load(exp))
        for fam in sorted(sel["oracle"].index.get_level_values(0).unique()):
            if exp in LFR_EXPS and (exp, fam) not in LFR_FAMS:
                continue
            r = dict(experiment=exp, family=fam)
            for prot, P in sel.items():
                # as the scorecard: each test graph seed averaged over the noise range, then
                # the ten seeds resampled (paired)
                g = P.loc[fam].dropna().groupby(level="seed").apply(
                    lambda h: pd.Series({m: np.trapezoid(h.sort_index(level="p")[m].to_numpy(),
                                                         h.sort_index(level="p").index.get_level_values("p").to_numpy(dtype=float))
                                         / (h.index.get_level_values("p").max() - h.index.get_level_values("p").min())
                                         for m in ("hc", "leiden", "leiden_j")}))
                for m in ("hc", "leiden", "leiden_j"):
                    r[f"{prot}_{m}"] = g[m].mean()
                diff = (g["hc"] - g["leiden_j"]).to_numpy()
                bs = [diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(2000)]
                r[f"{prot}_diff"], r[f"{prot}_lo"], r[f"{prot}_hi"] = diff.mean(), *np.percentile(bs, [2.5, 97.5])
                r[f"{prot}_n"] = len(diff)
            rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(DATA, "jaccard_control.csv"), index=False)

    def label(e, f):
        if e in LFR_EXPS:
            return LFR_FAMS[(e, f)], "", (2, list(LFR_FAMS).index((e, f)), "")
        return ring_name(e, f)

    T["lab"], T["on"], T["key"] = zip(*[label(e, f) for e, f in zip(T.experiment, T.family)])
    T = T.sort_values("key", key=lambda s: s.map(str)).reset_index(drop=True)
    T = T.iloc[sorted(range(len(T)), key=lambda i: T.key[i])]
    lines = [r"\begin{table*}[p]",
             r"\caption{Control for the commonality alone: Leiden-RB modularity on the graph weighted by "
             r"the closed-neighbourhood Jaccard commonality of each edge (Leiden-J, the 30 resolutions of "
             r"Leiden), against HyperCommon (HC) and unweighted Leiden. Mean $\Omega$ over the noise range "
             r"($p\in[0,0.8]$ for rings, $\mu\in[0.1,0.6]$ for LFR), ten test graphs per noise level; "
             r"oracle: best parameter per graph; label-free: largest overlapping modularity. HC $-$ L-J: "
             r"mean difference with the paired 95\% bootstrap interval over the ten test seeds, each averaged over the noise range. Weighted Leiden iterates until no improvement, at most 100 times. Bold: the "
             r"best of the three.}",
             r"\label{tab:S-jaccard-control}",
             r"\scriptsize\setlength{\tabcolsep}{3pt}",
             r"\begin{ruledtabular}",
             r"\begin{tabular}{llccc c|ccc c}",
             r" & & \multicolumn{4}{c|}{oracle} & \multicolumn{4}{c}{label-free} \\",
             r"Setting & $o_n$ & HC & Leiden & Leiden-J & HC $-$ L-J & HC & Leiden & Leiden-J & HC $-$ L-J \\",
             r"\hline"]
    prev = None
    for _, r in T.iterrows():
        cells = []
        for prot in ("oracle", "label_free"):
            v = [r[f"{prot}_{m}"] for m in ("hc", "leiden", "leiden_j")]
            b = max(v)
            cells += [(r"\textbf{%.2f}" if x >= b - 0.005 else "%.2f") % x for x in v]
            cells.append(f"${r[f'{prot}_diff']:+.2f}$ [${r[f'{prot}_lo']:+.2f},{r[f'{prot}_hi']:+.2f}$]")
        lab = r.lab if r.lab != prev else ""
        if r.key[0] == 2 and (prev is None or not str(prev).startswith("LFR")):
            lines.append(r"\hline")
        prev = r.lab
        lines.append(f"{lab} & {r.on} & " + " & ".join(cells) + r" \\")
    for prot in ("oracle", "label_free"):
        for kind, sel in (("rings", T.key.map(lambda k: k[0] < 2)), ("LFR", T.key.map(lambda k: k[0] == 2))):
            t = T[sel]
            print(f"{prot} {kind}: HC ahead of Leiden-J {(t[f'{prot}_lo'] > 0).sum()}, level "
                  f"{((t[f'{prot}_lo'] <= 0) & (t[f'{prot}_hi'] >= 0)).sum()}, behind {(t[f'{prot}_hi'] < 0).sum()} "
                  f"of {len(t)}; Leiden-J above Leiden in {(t[f'{prot}_leiden_j'] > t[f'{prot}_leiden'] + 0.005).sum()}")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    with open(TEX, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
