"""fig_omega_t.pdf: Omega vs threshold t for HyperCommon on overlapping uniform rings,
one curve per rewiring level p (mean over 10 test graphs), with the exact-recovery
window of Theorem 2 shaded for the noise-free case."""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from fig_real_networks import INK, INK2, GRID, FIG, panel_label

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "bench", "results", "rings_uniform.csv"))
df = df[(df.method == "hc") & (df.split == "test") & df.error.isna()]
df["t"] = df.param.astype(float)

def window(k):
    lo = 1 / (2 * k + 1)
    hi = min(k / (2 * k + 2), max((k - 2) / (k + 4), (k - 1) / (2 * k + 3)))
    return lo, hi

fams = [("k8_s50_o05", 8, r"$k=8,\ s=50,\ o_n=0.05$"), ("k16_s100_o05", 16, r"$k=16,\ s=100,\ o_n=0.05$")]
ps = [0.0, 0.2, 0.4, 0.6]
cmap = plt.get_cmap("viridis")
iso = pd.read_csv(os.path.join(HERE, "..", "bench", "results2", "rings_isolated.csv"))
iso = iso[(iso.method == "hc") & (iso.split == "test") & iso.error.isna() & (iso.p == 0)]
iso["t"] = iso.param.astype(float)
ISO = {8: "k8_s50_o01_iso", 16: "k16_s100_o01_iso"}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
for ax, (fam, k, title), letter in zip(axes, fams, "ab"):
    lo, hi = window(k)
    ax.axvspan(lo, hi, color="#eda100", alpha=0.18, lw=0)
    ax.text((lo + hi) / 2, 1.045, "predicted window, $p=0$ (isolated merges)", ha="center", fontsize=6.5, color=INK2)
    for i, p in enumerate(ps):
        g = df[(df.family == fam) & (df.p.round(2) == p)].groupby("t").omega.mean()
        c = cmap(0.85 * i / (len(ps) - 1))
        ax.plot(g.index, g.values, color=c, lw=1.6, label=f"$p={p:.1f}$")
    gi = iso[iso.family == ISO[k]].groupby("t").omega.mean()
    ax.plot(gi.index, gi.values, color=INK, lw=1.0, ls=(0, (3, 2)), label="isolated merges, $p=0$")
    ax.set_xlim(0, 0.6); ax.set_ylim(-0.02, 1.1)
    ax.set_xlabel("threshold $t$"); ax.set_title(title)
    ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True); panel_label(ax, letter)
axes[0].set_ylabel(r"$\Omega$ (mean of 10 graphs)")
axes[1].legend(frameon=False, fontsize=6.5, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig_omega_t.pdf")); fig.savefig(os.path.join(FIG, "fig_omega_t.png"), dpi=200)
