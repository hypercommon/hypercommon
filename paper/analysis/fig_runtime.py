"""
fig_runtime.pdf: running time against network size (bench/runtime.py output),
mean of three repetitions, log-log, with fitted power-law exponents for the
HyperCommon build and Leiden. Writes figures/data/runtime_summary.csv.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_real_networks import INK2, GRID, FIG, panel_label

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(HERE, "..", "bench", "results", "runtime.csv"))
g = d.groupby(["family", "n"]).mean(numeric_only=True).reset_index()
# the local query timed with its lazy node cache (bench/runtime_local.py) replaces the first timing
loc = os.path.join(HERE, "..", "bench", "results", "runtime_local.csv")
if os.path.exists(loc):
    L = pd.read_csv(loc).groupby(["family", "n"]).hc_local_s.mean()
    g = g.drop(columns="hc_local_s").merge(L.reset_index(), on=["family", "n"], how="left")
g.to_csv(os.path.join(HERE, "..", "figures", "data", "runtime_summary.csv"), index=False)

SERIES = [("hc_build_s", "HyperCommon: wedges + scores", "#2a78d6", "-"),
          ("hc_one_t_s", "HyperCommon: one threshold", "#2a78d6", "--"),
          ("hc_local_s", "HyperCommon: one local community", "#2a78d6", ":"),
          ("leiden_s", "Leiden", "#eb6834", "-"),
          ("infomap_s", "Infomap", "#e87ba4", "-")]
FAMS = [("rings_k8", r"rings, $k=8$"), ("rings_k16", r"rings, $k=16$"), ("lfr", r"LFR, $\langle k\rangle=20$")]

fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True)
fits = []
for ax, (fam, title), letter in zip(axes, FAMS, "abc"):
    h = g[g.family == fam].sort_values("n")
    for col, lab, c, ls in SERIES:
        ax.loglog(h.n, h[col], ls, color=c, lw=1.5, marker="o", ms=2.5, label=lab)
        if col in ("hc_build_s", "leiden_s", "infomap_s") and len(h) >= 3:
            sl = np.polyfit(np.log(h.n), np.log(h[col]), 1)[0]
            fits.append(dict(family=fam, series=col, exponent=sl))
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("nodes $n$")
    ax.grid(color=GRID, lw=0.5, which="both")
    panel_label(ax, letter)
axes[0].set_ylabel("seconds (one core)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=6.6, bbox_to_anchor=(0.5, 0.0), handlelength=4.5)
fig.tight_layout(rect=(0, 0.14, 1, 1))
fig.savefig(os.path.join(FIG, "fig_runtime.pdf"))
fig.savefig(os.path.join(FIG, "fig_runtime.png"), dpi=200)
pd.DataFrame(fits).to_csv(os.path.join(HERE, "..", "figures", "data", "runtime_exponents.csv"), index=False)
print(pd.DataFrame(fits).round(2).to_string(index=False))
