"""fig_real_family.pdf: commonality as a family on real networks.
Top: inner/outer box plots of four functions (each on its own scale).
Bottom: AUC of all eight functions, pooled (filled) and degree-stratified (open)."""
import json, os
import numpy as np, matplotlib.pyplot as plt
from fig_real_networks import mpl, C, INK, INK2, GRID, FIG, OUT, panel_label   # shared style

NETS = [("Karate", "Karate club"), ("Dblp", "DBLP"), ("Amazon", "Amazon"), ("Youtube", "YouTube")]
SHOW = [("jaccard", "J"), ("count", "I"), ("cosine", "cos"), ("min_norm", r"$I/\min$")]
ALL = [("jaccard", "J"), ("count", "I"), ("i2_over_u", r"$I^2/U$"),
       ("i_over_sqrtu", r"$I/U^{1/2}$ ($\equiv I^2/U$)"), ("cosine", "cos"), ("min_norm", r"$I/\min$"),
       ("max_norm", r"$I/\max$"), ("resource_alloc", "RA")]

fig = plt.figure(figsize=(7.0, 5.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1], hspace=0.45)
top = gs[0].subgridspec(1, 4, wspace=0.35)
for col, (key, title) in enumerate(NETS):
    d = json.load(open(os.path.join(OUT, f"{key}_family.json")))
    ax = fig.add_subplot(top[col])
    pos, data, cols = [], [], []
    for i, (f, lab) in enumerate(SHOW):
        b = d["f"][f]["box"]
        for j, g in enumerate(("inner", "outer")):
            q = b[g]
            # normalise each function by its inner 95th percentile so the four share an axis
            s = max(b["inner"][4], b["outer"][4]) or 1.0
            data.append([x / s for x in q]); pos.append(i * 3 + j); cols.append(C[g])
    for p_, q, c in zip(pos, data, cols):
        ax.plot([p_, p_], [q[0], q[1]], color=INK2, lw=0.7)
        ax.plot([p_, p_], [q[3], q[4]], color=INK2, lw=0.7)
        ax.add_patch(plt.Rectangle((p_ - 0.38, q[1]), 0.76, q[3] - q[1], color=c, alpha=0.85, lw=0))
        ax.plot([p_ - 0.38, p_ + 0.38], [q[2], q[2]], color=INK, lw=1.1)
    ax.set_xticks([i * 3 + 0.5 for i in range(len(SHOW))])
    ax.set_xticklabels([lab for _, lab in SHOW], rotation=35, fontsize=6.5)
    ax.set_xlim(-0.8, len(SHOW) * 3 - 1.2); ax.set_ylim(0, 1.05)
    ax.set_title(title); ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
    if col == 0:
        ax.set_ylabel("f on edges\n(scaled to its 95th pct.)")
    panel_label(ax, "abcd"[col])
from matplotlib.patches import Patch
fig.legend([Patch(color=C["inner"]), Patch(color=C["outer"])], ["inner (shared recorded group)", "outer (labelled, no shared group)"],
           loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))

ax = fig.add_subplot(gs[1])
x = np.arange(len(ALL))
for k, (key, title) in enumerate(NETS):
    d = json.load(open(os.path.join(OUT, f"{key}_family.json")))
    off = (k - 1.5) * 0.18
    pooled = [d["f"][f]["auc"] for f, _ in ALL]
    strat = [d["f"][f]["auc_strat_mean"] for f, _ in ALL]
    color = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"][k]
    mk = "osD^"[k]
    ax.plot(x + off, pooled, mk, color=color, ms=4.5, label=title)
    if strat[0] is not None:
        ax.plot(x + off, strat, mk, mfc="white", mec=color, ms=4.5, mew=1.1)
ax.axhline(0.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
ax.set_xticks(x); ax.set_xticklabels([lab for _, lab in ALL], fontsize=7)
ax.set_ylabel("AUC, inner vs\nouter edges"); ax.set_ylim(0.42, 0.92)
ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
ax.plot([], [], "o", color=INK2, ms=4.5, label="all edges"); ax.plot([], [], "o", mfc="white", mec=INK2, ms=4.5, label="within degree strata")
ax.legend(ncol=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), fontsize=6.8)
panel_label(ax, "e")
fig.savefig(os.path.join(FIG, "fig_real_family.pdf")); fig.savefig(os.path.join(FIG, "fig_real_family.png"), dpi=200)
