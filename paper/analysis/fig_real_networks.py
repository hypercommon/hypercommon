"""
Figures for the empirical-commonality section.

fig_real_commonality.pdf  (two rows, four datasets)
  (a-d) closed Jaccard of edges by class: inner / outer / unlabelled
  (e-h) AUC(inner vs outer) within max-degree strata, with the pooled AUC and
        the AUC the thesis effectively computed (inner vs outer+unlabelled)

fig_real_triples.pdf  precision-recall of "all three nodes share a community"
  for the pair rule (both wedge edges >= t) vs the triple rule (all three
  pairs >= t), with the closed-triangle-only point.
"""

import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

DATASETS = [("Karate", "Karate club"), ("Dblp", "DBLP"), ("Amazon", "Amazon"), ("Youtube", "YouTube")]
C = {"inner": "#2a78d6", "outer": "#eb6834", "unlabelled": "#1baf7a"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "Times New Roman"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def load(name):
    with open(os.path.join(OUT, f"{name}.json")) as fh:
        a = json.load(fh)
    with open(os.path.join(OUT, f"{name}_v2.json")) as fh:
        b = json.load(fh)
    with open(os.path.join(OUT, f"{name}_edge_sample.json")) as fh:
        s = json.load(fh)
    return a, b, s


def panel_label(ax, s):
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def fig_commonality():
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.9), gridspec_kw={"height_ratios": [1.15, 1]})
    letters = "abcdefgh"
    for col, (key, title) in enumerate(DATASETS):
        a, b, s = load(key)
        j = np.array(s["j"])
        cls = np.array(s["cls"])
        ax = axes[0, col]
        groups = [g for g in ("inner", "outer", "unlabelled") if (cls == g).sum() > 0]
        data = [j[cls == g] for g in groups]
        bp = ax.boxplot(data, positions=range(len(groups)), widths=0.55, showfliers=False,
                        patch_artist=True, medianprops={"color": INK, "lw": 1.1},
                        whiskerprops={"color": INK2, "lw": 0.7}, capprops={"color": INK2, "lw": 0.7},
                        boxprops={"lw": 0.0})
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(C[g])
            patch.set_alpha(0.85)
        counts = a["edges"]["class_counts"]
        ax.set_xticks(range(len(groups)))
        # thin spaces as digit separators (935 563, not 935,563)
        ax.set_xticklabels([f"{g}\n" + "$" + f"{counts[g]:,}".replace(",", r"\,") + "$" for g in groups],
                           fontsize=6.3)
        ax.set_title(title, color=INK)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        if col == 0:
            ax.set_ylabel(r"closed Jaccard $J(u,v)$, edges")
        panel_label(ax, letters[col])

        ax = axes[1, col]
        strata = a["edges"]["degree_stratified"]
        ax.axhline(0.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
        if strata:
            x = np.arange(len(strata))
            ax.plot(x, [st["auc"] for st in strata], "-o", color=C["inner"], lw=1.6, ms=3.5,
                    label="within degree stratum")
            ax.set_xticks(x)
            ax.set_xticklabels([f"{st['dmax_lo']}" + (f"–{st['dmax_hi']}" if st["dmax_hi"] < 1e8 else "+")
                                for st in strata], rotation=40, fontsize=6.3)
            xr = (-0.4, len(strata) - 0.6)
        else:
            xr = (-0.5, 0.5)
            ax.set_xticks([])
            ax.text(0, 0.64, "too few outer edges\n(11) to stratify", ha="center", fontsize=6.5, color=INK2)
        ax.axhline(a["edges"]["auc_inner_vs_outer"], color=C["inner"], lw=0.9, ls="--", label="pooled, labelled")
        if counts["unlabelled"]:
            ax.axhline(a["edges"]["auc_inner_vs_pooled"], color=C["unlabelled"], lw=0.9, ls="-.",
                       label="pooled, unlabelled as outer")
        ax.set_xlim(*xr)
        ax.set_ylim(0.35, 0.95)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlabel(r"larger endpoint degree $d_{\max}$", fontsize=7)
        if col == 0:
            ax.set_ylabel("AUC (inner vs outer)")
        panel_label(ax, letters[4 + col])
    h, l = axes[1, 1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(h_pad=1.2, w_pad=0.8, rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(FIG, "fig_real_commonality.pdf"))
    fig.savefig(os.path.join(FIG, "fig_real_commonality.png"), dpi=200)


def fig_triples():
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 3.6), gridspec_kw={"height_ratios": [1.5, 1]})
    for col, (key, title) in enumerate(DATASETS[1:]):
        _, b, _ = load(key)
        mr = b["wedges"]["matched_recall"]
        ax = axs[1, col]
        R = sorted(mr, key=float)
        x = [float(r) for r in R]
        d = np.array([mr[r]["diff"] for r in R])
        lo = np.array([mr[r]["diff_lo"] for r in R])
        hi = np.array([mr[r]["diff_hi"] for r in R])
        ax.axhline(0, color=INK2, lw=0.6, ls=(0, (2, 2)))
        ax.errorbar(x, d, yerr=[d - lo, hi - d], fmt="o-", color=C["inner"], ms=3, lw=1.2, capsize=2)
        ax.set_xlim(0, 1)
        ax.set_xlabel("recall (pure triples admitted)")
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        panel_label(ax, "def"[col])
    axs[1, 0].set_ylabel("precision difference\n(triple $-$ pair)")
    axes = axs[0]
    for ax, (key, title), letter in zip(axes, DATASETS[1:], "abc"):
        _, b, _ = load(key)
        w = b["wedges"]
        for rule, color, lab in (("pair", C["outer"], r"pair rule: every edge $\geq t$"),
                                 ("triple", C["inner"], r"triple rule: all three pairs $\geq t$")):
            pts = [(p["recall"], p["precision"]) for p in w[rule] if p["precision"] is not None and p["kept"] >= 200]
            r, p = zip(*pts)
            ax.plot(r, p, color=color, lw=1.6, label=lab)
        tri = w["triangle_only"]
        ax.plot(tri["recall"], tri["precision"], "D", color=INK, ms=4, label="closed triangles only")
        ax.axhline(w["base_rate"], color=INK2, lw=0.6, ls=(0, (2, 2)), label="base rate")
        ax.set_xlim(0, 1)
        ax.set_title(title)
        ax.set_xlabel("recall (pure triples admitted)")
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        panel_label(ax, letter)
    axes[0].set_ylabel("precision (triple within\none recorded group)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(w_pad=1.0, h_pad=0.8, rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(FIG, "fig_real_triples.pdf"))
    fig.savefig(os.path.join(FIG, "fig_real_triples.png"), dpi=200)


if __name__ == "__main__":
    import sys
    if "triples" not in sys.argv[1:]:
        fig_commonality()
    fig_triples()
