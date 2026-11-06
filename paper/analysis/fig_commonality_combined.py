"""fig_commonality.pdf: the main-text figure of the empirical-commonality section.

Reads only the cached results in analysis/out (no recomputation on the SNAP graphs).

(a) AUC(inner vs outer) within strata of the larger endpoint degree d_max, for DBLP,
    Amazon and YouTube (points), with the pooled AUC (dashed) and the AUC with
    unlabelled edges counted as outer (dash-dotted)       [old Fig. 1(f-h)]
(b) AUC of the eight commonality functions, all labelled edges (filled) and within
    strata of both endpoint degrees (open), four networks [old Fig. 2(e)]
(c) precision difference triple rule - pair rule at equal recall, with 95% intervals
                                                          [old Fig. 3(d-f)]
The full versions of the three old figures stay in the Supplemental Material.
"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from fig_real_networks import mpl, INK, INK2, GRID, FIG, OUT, panel_label  # shared style

NETS = [("Karate", "Karate club"), ("Dblp", "DBLP"), ("Amazon", "Amazon"), ("Youtube", "YouTube")]
COL = {"Karate": "#2a78d6", "Dblp": "#eb6834", "Amazon": "#1baf7a", "Youtube": "#6250d6"}
MK = {"Karate": "o", "Dblp": "s", "Amazon": "D", "Youtube": "^"}
ALL = [("jaccard", "$J$"), ("count", "$I$"), ("i2_over_u", r"$I^2/U$"),
       ("i_over_sqrtu", r"$I/U^{1/2}$"), ("cosine", "cos"), ("min_norm", r"$I/\min$"),
       ("max_norm", r"$I/\max$"), ("resource_alloc", "RA")]


def jload(name):
    with open(os.path.join(OUT, name)) as fh:
        return json.load(fh)


fig = plt.figure(figsize=(7.0, 2.55))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.45, 1.0], wspace=0.34)

# (a) stratified AUC against d_max
ax = fig.add_subplot(gs[0])
ax.axhline(0.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
labels = None
for k, (key, title) in enumerate(NETS[1:]):
    e = jload(f"{key}.json")["edges"]
    strata = e["degree_stratified"]
    x = np.arange(len(strata)) + (k - 1) * 0.12
    ax.plot(x, [st["auc"] for st in strata], "-", marker=MK[key], color=COL[key], lw=1.2, ms=3.5,
            label=title)
    ax.axhline(e["auc_inner_vs_outer"], color=COL[key], lw=0.8, ls="--")
    ax.axhline(e["auc_inner_vs_pooled"], color=COL[key], lw=0.8, ls="-.")
    labels = [f"{st['dmax_lo']}" + (f"–{st['dmax_hi']}" if st["dmax_hi"] < 1e8 else "+") for st in strata]
ax.set_xticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=40, fontsize=6.5)
ax.set_xlim(-0.4, len(labels) - 0.6)
ax.set_ylim(0.38, 0.92)
ax.set_xlabel(r"larger endpoint degree $d_{\max}$")
ax.set_ylabel("AUC, inner vs outer edges")
ax.grid(axis="y", color=GRID, lw=0.5)
ax.set_axisbelow(True)
panel_label(ax, "a")

# (b) the family of functions
ax = fig.add_subplot(gs[1])
x = np.arange(len(ALL))
for k, (key, title) in enumerate(NETS):
    d = jload(f"{key}_family.json")
    off = (k - 1.5) * 0.18
    pooled = [d["f"][f]["auc"] for f, _ in ALL]
    strat = [d["f"][f]["auc_strat_mean"] for f, _ in ALL]
    ax.plot(x + off, pooled, MK[key], color=COL[key], ms=3.8, label=title)
    if strat[0] is not None:
        ax.plot(x + off, strat, MK[key], mfc="white", mec=COL[key], ms=3.8, mew=1.0)
ax.axhline(0.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
ax.set_xticks(x)
ax.set_xticklabels([lab for _, lab in ALL], fontsize=6.8, rotation=30)
ax.set_xlim(-0.5, len(ALL) - 0.5)
ax.set_ylim(0.42, 0.95)
ax.set_xlabel("commonality function")
ax.grid(axis="y", color=GRID, lw=0.5)
ax.set_axisbelow(True)
panel_label(ax, "b")

# (c) triple rule minus pair rule at equal recall
ax = fig.add_subplot(gs[2])
ax.axhline(0, color=INK2, lw=0.6, ls=(0, (2, 2)))
for k, (key, title) in enumerate(NETS[1:]):
    mr = jload(f"{key}_v2.json")["wedges"]["matched_recall"]
    R = sorted(mr, key=float)
    xr = np.array([float(r) for r in R]) + (k - 1) * 0.008
    dd = np.array([mr[r]["diff"] for r in R])
    lo = np.array([mr[r]["diff_lo"] for r in R])
    hi = np.array([mr[r]["diff_hi"] for r in R])
    ax.errorbar(xr, dd, yerr=[dd - lo, hi - dd], fmt="-", marker=MK[key], color=COL[key], ms=3.5,
                lw=1.1, capsize=1.8, elinewidth=0.8, label=title)
ax.set_xlim(0.05, 0.45)
ax.set_xticks([0.1, 0.2, 0.3, 0.4])
ax.set_xlabel("recall (pure triples admitted)")
ax.set_ylabel("precision, triple $-$ pair rule")
ax.grid(color=GRID, lw=0.5)
ax.set_axisbelow(True)
panel_label(ax, "c")

# one legend for the whole figure: networks (same colour and marker in every panel), then line keys
h = [Line2D([], [], color=COL[k], marker=MK[k], ms=3.8, lw=1.1 if k != "Karate" else 0) for k, _ in NETS]
l = [t for _, t in NETS]
h += [Line2D([], [], color=INK2, lw=0.8, ls="--"), Line2D([], [], color=INK2, lw=0.8, ls="-."),
      Line2D([], [], color=INK2, marker="o", ms=3.8, lw=0),
      Line2D([], [], color=INK2, marker="o", mfc="white", ms=3.8, lw=0)]
l += ["(a) pooled", "(a) unlabelled counted as outer", "(b) all labelled edges", "(b) within degree strata"]
fig.legend(h, l, loc="upper center", ncol=8, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=6.6,
           handlelength=1.8, columnspacing=1.0, handletextpad=0.4)
fig.subplots_adjust(top=0.82)

fig.savefig(os.path.join(FIG, "fig_commonality.pdf"))
fig.savefig(os.path.join(FIG, "fig_commonality.png"), dpi=200)
