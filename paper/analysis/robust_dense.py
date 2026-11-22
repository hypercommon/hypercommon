"""
Robustness to noise on the ring lattices, p = 0, 0.01, ..., 1: the main grid (results2/<exp>.csv,
p in steps of 0.05) plus the dense sweep (bench/run2.py, run_oslom.py and run_sbm.py with
DENSE_P=1, results2/<exp>_dense*.csv), test graphs only, every method of the main runs
including OSLOM and the stochastic block model (ablations and native-default rows excluded).

Per family, method and protocol (oracle: best parameter per graph; label-free: largest
overlapping modularity), the mean Omega over the ten test graphs at each p with a 95%
bootstrap band, and from that curve
  p_drop(L)  the first p at which the mean falls below L (linear interpolation between grid
             points; inf if it never does), for L = 0.5 ... 0.9
  area       the area under the curve over p in [0, 1] (the mean Omega over the noise range)
Planted structure left at each p (bench/robust_structure.py):
  p_weak     the first p at which fewer than half of the planted rings are communities in the
             weak sense (more edge ends inside than outside)
  p_clust    the first p at which the clustering coefficient falls below 5% of its p = 0 value

    python robust_dense.py
Writes figures/data/robust_curves.csv, robust_drop.csv, robust_structure.csv and
figures/fig_robust.(pdf|png) (six families) and fig_robust_all.(pdf|png) (all 28).
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_real_networks import GRID, FIG, panel_label

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.join(HERE, "..", "bench", "results2")
DATA = os.path.join(HERE, "..", "figures", "data")
EXPS = ["rings_uniform", "rings_hetero", "rings_hetero_overlap", "rings_pref"]
DROP = {"hc_2test", "hc_tri", "hc_econn", "leiden_g1", "leiden_cpm_dens", "linkcomm_pd"}
LEVELS = (0.5, 0.6, 0.7, 0.8, 0.9)
LEVEL_FIG = 0.6
PANELS = [  # the five ring panels of fig_benchmarks, and the equal-size unequal-degree family
    ("rings_uniform", "k8_s50_o10", r"uniform rings, $k=8$, $o_n=0.10$"),
    ("rings_uniform", "k16_s100_o05", r"uniform rings, $k=16$, $o_n=0.05$"),
    ("rings_hetero", "two_worlds", "two worlds"),
    ("rings_hetero", "three_big_many_small", "three big, many small"),
    ("rings_hetero", "gradient10", "gradient of ten rings"),
    ("rings_hetero", "k_only_s100", "equal size, unequal degree"),
]
COLOR = {"hc": "#0072B2", "linkcomm": "#7a7f87", "cpm": "#009E73", "leiden": "#D55E00",
         "leiden_cpm": "#CC79A7", "infomap": "#56B4E9", "slpa": "#E69F00", "oslom": "#000000",
         "sbm": "#882255", "walktrap": "#999933", "demon": "#44AA99", "angel": "#DDCC77",
         "lpa": "#AA4499", "scan": "#88CCEE", "scan_j": "#117733",
         "bigclam": "#661100", "egosplit": "#6699CC"}
LABEL = {"hc": "HC", "linkcomm": "Link", "cpm": "CPM", "leiden": "Leiden", "leiden_cpm": "L-CPM",
         "infomap": "Infomap", "slpa": "SLPA", "oslom": "OSLOM", "sbm": "SBM", "walktrap": "Walktrap",
         "demon": "DEMON", "angel": "ANGEL", "lpa": "LPA", "scan": "SCAN", "scan_j": "SCAN-J",
         "bigclam": "BigClam", "egosplit": "Ego"}
rng = np.random.default_rng(11)


def load(exp):
    parts = []
    for suf in ("", "_dense"):
        for tag in ("", "_oslom", "_sbm"):
            f = os.path.join(R2, f"{exp}{suf}{tag}.csv")
            if os.path.exists(f) and os.path.getsize(f) > 0:
                parts.append(pd.read_csv(f, usecols=["family", "p", "seed", "split", "method", "omega", "eq"],
                                         low_memory=False))
    d = pd.concat(parts, ignore_index=True)
    d = d[(d.split == "test") & d.omega.notna() & ~d.method.isin(DROP)].copy()
    d["p"] = d.p.astype(float).round(2)
    return d


def curves():
    out = []
    for exp in EXPS:
        d = load(exp)
        for prot, col in (("oracle", "omega"), ("label_free", "eq")):
            x = d[d[col].notna()]
            sel = x.loc[x.groupby(["family", "method", "p", "seed"])[col].idxmax()]
            for (fam, m, p), g in sel.groupby(["family", "method", "p"]):
                v = g.omega.to_numpy()
                b = v[rng.integers(0, len(v), (1000, len(v)))].mean(1)
                out.append(dict(experiment=exp, family=fam, method=m, protocol=prot, p=p, omega=v.mean(),
                                lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5), graphs=len(v)))
    return pd.DataFrame(out)


def first_below(ps, ys, level):
    for i in range(len(ps)):
        if ys[i] < level:
            if i == 0:
                return ps[0]
            a, b = ys[i - 1], ys[i]
            return ps[i - 1] + (a - level) / (a - b) * (ps[i] - ps[i - 1])
    return np.inf


def structure():
    f = os.path.join(R2, "robust_structure.csv")
    if not os.path.exists(f):
        return pd.DataFrame(columns=["experiment", "family", "p_weak", "p_clust"]), None
    s = pd.read_csv(f)
    s["p"] = s.p.round(2)
    m = s.groupby(["experiment", "family", "p"])[["inside", "weak", "strong", "clustering"]].mean().reset_index()
    rows = []
    for (e, fam), g in m.groupby(["experiment", "family"]):
        g = g.sort_values("p")
        c0 = g.clustering.iloc[0]
        rows.append(dict(experiment=e, family=fam,
                         p_weak=first_below(g.p.values, g.weak.values, 0.5),
                         p_inside=first_below(g.p.values, g.inside.values, 0.5),
                         p_clust=first_below(g.p.values, g.clustering.values / c0, 0.05)))
    return pd.DataFrame(rows), m


def drops(C):
    rows = []
    for (e, fam, prot, m), g in C.groupby(["experiment", "family", "protocol", "method"]):
        g = g.sort_values("p")
        ps, ys = g.p.values, g.omega.values
        r = dict(experiment=e, family=fam, protocol=prot, method=m, p_max=ps.max(), n_p=len(ps),
                 area=float(np.trapezoid(ys, ps) / (ps.max() - ps.min())))
        r.update({f"p_drop_{L}": first_below(ps, ys, L) for L in LEVELS})
        rows.append(r)
    D = pd.DataFrame(rows)
    # HyperCommon against the best other method per family and protocol, by each score
    comp = []
    for (e, fam, prot), g in D.groupby(["experiment", "family", "protocol"]):
        g = g.set_index("method")
        if "hc" not in g.index:
            continue
        # compare only methods whose curve reaches the family's full range (a method still
        # running would be scored on a shorter range, which favours it)
        full = g.p_max.max()
        late = sorted(g.index[g.p_max < full - 1e-9])
        g = g[g.p_max >= full - 1e-9]
        o = g.drop(index="hc")
        r = dict(experiment=e, family=fam, protocol=prot, p_max=full, left_out=",".join(late),
                 hc_area=g.loc["hc", "area"], other_area=o.area.max(), other_area_by=o.area.idxmax())
        for L in LEVELS:
            k = f"p_drop_{L}"
            r[f"hc_{k}"] = g.loc["hc", k]
            r[f"other_{k}"] = o[k].max()
            r[f"other_{k}_by"] = o[k].idxmax()
            r[f"hc_rank_{k}"] = int((o[k] > g.loc["hc", k] + 1e-9).sum()) + 1
        comp.append(r)
    return D, pd.DataFrame(comp)


def draw(ax, C, S, exp, fam, title, letter=None, small=False):
    g = C[(C.experiment == exp) & (C.family == fam)]
    if g.empty:
        return
    area = g.groupby("method").omega.mean()
    others = area.drop(index="hc", errors="ignore").sort_values(ascending=False).index[:3]
    st = S[(S.experiment == exp) & (S.family == fam)] if len(S) else S
    if len(st) and np.isfinite(st.p_weak.iloc[0]):
        ax.axvspan(st.p_weak.iloc[0], 1.0, color="#bbbbbb", alpha=0.25, lw=0)
    if len(st) and np.isfinite(st.p_clust.iloc[0]):
        ax.axvline(st.p_clust.iloc[0], color="#555555", lw=0.8, ls=":")
    ax.axhline(LEVEL_FIG, color="#888888", lw=0.6, ls="--")
    for m in list(others)[::-1] + ["hc"]:
        h = g[g.method == m].sort_values("p")
        if h.empty:
            continue
        ax.fill_between(h.p, h.lo, h.hi, color=COLOR.get(m, "#444444"), alpha=0.15, lw=0)
        ax.plot(h.p, h.omega, color=COLOR.get(m, "#444444"), lw=2.0 if m == "hc" else 1.1, label=LABEL.get(m, m))
    ax.set_title(title, fontsize=6.5 if small else 7.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=5 if small else 5.8, frameon=False, loc="lower left", handlelength=1.4)
    if letter:
        panel_label(ax, letter)


def figure(C, S, protocol):
    C = C[C.protocol == protocol]
    ylab = r"$\Omega$ (label-free parameter)" if protocol == "label_free" else r"$\Omega$ (best parameter)"
    suf = "" if protocol == "label_free" else "_oracle"
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.5), sharey=True)
    for ax, (exp, fam, title), letter in zip(axes.ravel(), PANELS, "abcdef"):
        draw(ax, C, S, exp, fam, title, letter)
        ax.set_xlabel("rewired fraction $p$", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylab)
    fig.tight_layout(h_pad=1.0)
    fig.savefig(os.path.join(FIG, f"fig_robust{suf}.pdf"))
    fig.savefig(os.path.join(FIG, f"fig_robust{suf}.png"), dpi=200)
    plt.close(fig)
    fams = C[["experiment", "family"]].drop_duplicates().values.tolist()
    nc = 4
    nr = int(np.ceil(len(fams) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(7.0, 1.55 * nr), sharex=True, sharey=True)
    for ax, (exp, fam) in zip(axes.ravel(), fams):
        draw(ax, C, S, exp, fam, f"{exp.replace('rings_', '')}: {fam}", small=True)
    for ax in axes.ravel()[len(fams):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("$p$", fontsize=7)
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    fig.savefig(os.path.join(FIG, f"fig_robust_all{suf}.pdf"))
    fig.savefig(os.path.join(FIG, f"fig_robust_all{suf}.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    C = curves()
    C.to_csv(os.path.join(DATA, "robust_curves.csv"), index=False)
    S, Sm = structure()
    if Sm is not None:
        Sm.to_csv(os.path.join(DATA, "robust_structure.csv"), index=False)
    D, comp = drops(C)
    D.to_csv(os.path.join(DATA, "robust_drop_by_method.csv"), index=False)
    comp = comp.merge(S, on=["experiment", "family"], how="left")
    comp.to_csv(os.path.join(DATA, "robust_drop.csv"), index=False)
    for prot in ("label_free", "oracle"):
        figure(C, S, prot)
    k = f"p_drop_{LEVEL_FIG}"
    cols = ["experiment", "family", "p_max", "left_out", f"hc_{k}", f"other_{k}", f"other_{k}_by", f"hc_rank_{k}",
            "hc_area", "other_area", "other_area_by", "p_weak", "p_clust"]
    with pd.option_context("display.width", 250):
        for prot in ("label_free", "oracle"):
            print(prot)
            print(comp[comp.protocol == prot][[c for c in cols if c in comp]].round(3).to_string(index=False))
