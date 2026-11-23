"""
The robustness table (Table S-robust): for each of the 28 ring settings and both protocols,
HyperCommon's p_0.6 (the first p at which its mean Omega over the ten test graphs falls below
0.6), the other method that holds longest and its p_0.6, HyperCommon's rank, and HyperCommon
minus that method; the paired 95% bootstrap interval is of HC minus the best of ALL other methods,
every mean curve recomputed on the same resample of the ten test graphs and the best other method
chosen again in each (a curve that never falls below 0.6 counts as 1). The same for the area under
the mean curve (area_diff, area_lo, area_hi in the csv).
Also the planted-structure markers p_weak and p_clust. Reads what robust_dense.py writes, and
the per-graph rows through its loader.

    python robust_table.py
Writes figures/data/robust_table.csv and sections/table_robust_dense.tex.
"""
import os

import numpy as np
import pandas as pd

from robust_dense import EXPS, load, first_below, DATA, LABEL

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, "..", "sections", "table_robust_dense.tex")
L = 0.6
rng = np.random.default_rng(5)

FAM = {  # family -> (name, o_n); the overlap level is in the name for the overlap sweeps
    "k8_s50": ("uniform, $k=8$", None), "k16_s100": ("uniform, $k=16$", None),
    "dense_even": ("equal $k$, unequal $s$, $k=16$", "0.05"), "gradient10": ("gradient", "0.05"),
    "k_only_s100": ("equal $s$, unequal $k$", "0.05"), "sizes_only_k8": ("equal $k$, unequal $s$, $k=8$", "0.05"),
    "three_big_many_small": ("three big, many small", "0.05"), "two_worlds": ("two worlds", "0.05"),
}
ORDER = ["k8_s50", "k16_s100", "dense_even", "k_only_s100", "sizes_only_k8", "gradient10",
         "three_big_many_small", "two_worlds"]


def name(exp, fam):
    """(name, o_n, rewiring) of a setting."""
    if exp == "rings_pref":
        base = "uniform, $k=8$" if fam.startswith("k8") else "two worlds"
        a = "2" if fam.endswith("2") else "1"
        return base, "0.05", rf"$\alpha={a}$"
    for key, (nm, on) in FAM.items():
        if fam.startswith(key):
            rest = fam[len(key):]
            if rest.startswith("_o"):
                on = {"00": "0", "05": "0.05", "10": "0.10"}[rest[2:4]]
            return nm, on, "uniform"
    return fam, "", ""


def key(exp, fam):
    nm, on, rw = name(exp, fam)
    base = next(k for k in ORDER if fam.startswith(k))
    return (rw != "uniform", ORDER.index(base), on, rw)


def per_graph(exp):
    """Per test graph, method, protocol and p: Omega of the selected cover."""
    d = load(exp)
    out = []
    for prot, col in (("oracle", "omega"), ("label_free", "eq")):
        x = d[d[col].notna()]
        sel = x.loc[x.groupby(["family", "method", "p", "seed"])[col].idxmax()]
        out.append(sel[["family", "method", "p", "seed", "omega"]].assign(experiment=exp, protocol=prot))
    return pd.concat(out, ignore_index=True)


def boot_diff(g, a, b, n=2000):
    """HC minus method b in p_0.6 of the mean curve, resampling the test graphs (paired)."""
    W = {m: g[g.method == m].pivot_table(index="seed", columns="p", values="omega") for m in (a, b)}
    seeds = W[a].index.intersection(W[b].index).to_numpy()
    A, B = W[a].loc[seeds], W[b].loc[seeds]
    pa, pb = A.columns.to_numpy(), B.columns.to_numpy()
    va, vb = A.to_numpy(), B.to_numpy()
    out = np.empty(n)
    for i in range(n):
        idx = rng.integers(0, len(seeds), len(seeds))
        xa = min(first_below(pa, np.nanmean(va[idx], 0), L), 1.0)
        xb = min(first_below(pb, np.nanmean(vb[idx], 0), L), 1.0)
        out[i] = xa - xb
    return np.percentile(out, [2.5, 97.5])


def boot_all(g, others, n=2000):
    """HC minus the BEST of all other methods, recomputed in every resample of the test graphs
    (paired: the same seeds for every method), for p_0.6 and for the area under the mean curve.
    Returns (lo, hi) of each, and the share of resamples in which HC leads by each."""
    W = {m: g[g.method == m].pivot_table(index="seed", columns="p", values="omega")
         for m in ["hc"] + list(others)}
    seeds = W["hc"].index
    for m in others:
        seeds = seeds.intersection(W[m].index)
    seeds = seeds.to_numpy()
    P = {m: W[m].columns.to_numpy(dtype=float) for m in W}
    V = {m: W[m].loc[seeds].to_numpy() for m in W}

    def stats(idx):
        out = {}
        for m in W:
            y = np.nanmean(V[m][idx], 0)
            ps = P[m]
            out[m] = (min(first_below(ps, y, L), 1.0), float(np.trapezoid(y, ps) / (ps.max() - ps.min())))
        return out

    dp, da = np.empty(n), np.empty(n)
    for i in range(n):
        st = stats(rng.integers(0, len(seeds), len(seeds)))
        dp[i] = st["hc"][0] - max(st[m][0] for m in others)
        da[i] = st["hc"][1] - max(st[m][1] for m in others)
    full = stats(np.arange(len(seeds)))
    return dict(lo_all=np.percentile(dp, 2.5), hi_all=np.percentile(dp, 97.5),
                area_diff=full["hc"][1] - max(full[m][1] for m in others),
                area_lo=np.percentile(da, 2.5), area_hi=np.percentile(da, 97.5))


def main():
    comp = pd.read_csv(os.path.join(DATA, "robust_drop.csv"))
    G = pd.concat([per_graph(e) for e in EXPS], ignore_index=True)
    k = f"p_drop_{L}"
    rows = []
    for _, r in comp.iterrows():
        g = G[(G.experiment == r.experiment) & (G.family == r.family) & (G.protocol == r.protocol)]
        bo = r[f"other_{k}_by"]
        lo, hi = boot_diff(g, "hc", bo)
        late = set(str(r.get("left_out", "") or "").split(",")) - {"", "nan"}
        others = sorted(set(g.method) - {"hc"} - late)
        ba = boot_all(g, others)
        rows.append(dict(experiment=r.experiment, family=r.family, protocol=r.protocol,
                         hc=r[f"hc_{k}"], best_other=bo, best_other_p=r[f"other_{k}"],
                         rank=int(r[f"hc_rank_{k}"]), diff=min(r[f"hc_{k}"], 1) - min(r[f"other_{k}"], 1), lo=lo, hi=hi,
                         hc_area=r.hc_area, other_area=r.other_area, other_area_by=r.other_area_by,
                         p_weak=r.p_weak, p_clust=r.p_clust, left_out=r.get("left_out", ""), **ba))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(DATA, "robust_table.csv"), index=False)

    def f(x):
        return r"$>1$" if not np.isfinite(x) else f"{x:.2f}"

    lines = [
        r"\begin{table*}[p]",
        r"\caption{Robustness to rewiring on the 28 ring settings: $p_{0.6}$, the rewired "
        r"fraction at which the mean $\Omega$ over the ten test graphs first falls below $0.6$ "
        r"(linearly interpolated between grid points), for HyperCommon (HC) and for the other "
        r"method that holds longest, HC's rank among the seventeen methods and HC minus the best "
        r"other method, with the paired 95\% bootstrap interval of HC minus the best of all other "
        r"methods (the test graphs resampled, every mean curve recomputed and the best other "
        r"method chosen again in each resample). Label-free: every "
        r"parameter by the largest overlapping modularity; oracle: the best parameter per "
        r"graph. $p_{\rm weak}$, $p_{\rm clust}$: the planted-structure markers of "
        r"Sec.~\ref{sec:robust}; $p_{\rm clust}>1$: the clustering never falls below $5\%$ of "
        r"its value at $p=0$. Rewiring is uniform except where marked preferential (pref.) "
        r"with exponent $\alpha$.}",
        r"\label{tab:S-robust}",
        r"\scriptsize\setlength{\tabcolsep}{2pt}",
        r"\begin{ruledtabular}",
        r"\begin{tabular}{llcc|cl r c|cl r c}",
        r" & & & & \multicolumn{4}{c|}{label-free} & \multicolumn{4}{c}{oracle} \\",
        r"Setting & $o_n$ & $p_{\rm weak}$ & $p_{\rm clust}$ & HC & best other & "
        r"rank & HC $-$ best & HC & best other & rank & HC $-$ best \\",
        r"\hline",
    ]
    settings = sorted({(e, fm) for e, fm in zip(T.experiment, T.family)}, key=lambda s: key(*s))
    for e, fm in settings:
        nm, on, rw = name(e, fm)
        q = T[(T.experiment == e) & (T.family == fm)].set_index("protocol")
        cells = []
        for prot in ("label_free", "oracle"):
            x = q.loc[prot]
            cells.append(f"{f(x.hc)} & {LABEL[x.best_other]} {f(x.best_other_p)} & {x['rank']} & "
                         f"${x['diff']:+.2f}$ [${x.lo_all:+.2f},{x.hi_all:+.2f}$]")
        x = q.iloc[0]
        nm = nm if rw == "uniform" else f"{nm}, pref. {rw}"
        lines.append(f"{nm} & {on} & {f(x.p_weak)} & {f(x.p_clust)} & "
                     f"{cells[0]} & {cells[1]} \\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    with open(TEX, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    for prot in ("label_free", "oracle"):
        t = T[T.protocol == prot]
        sig_win = ((t.lo_all > 0)).sum()
        sig_loss = ((t.hi_all < 0)).sum()
        print(f"{prot}: rank 1 in {(t['rank'] == 1).sum()}, top 3 in {(t['rank'] <= 3).sum()} of {len(t)}; "
              f"HC ahead of the best other with the interval above 0 in {sig_win}, behind with it below 0 "
              f"in {sig_loss}, level in {len(t) - sig_win - sig_loss} (against the best of all others); "
              f"by area ahead {(t.area_lo > 0).sum()}, behind {(t.area_hi < 0).sum()}")
    with pd.option_context("display.width", 250):
        print(T.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
