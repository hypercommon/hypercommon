"""
fig_benchmarks.pdf: Omega vs noise for representative graph families, each method
at its per-graph best parameter (oracle protocol), mean over 10 test graphs with
95% bootstrap bands. Also writes figures/data/benchmarks_curves.csv.
"""
import json
import os

import matplotlib.pyplot as plt
import pandas as pd

from fig_real_networks import INK2, GRID, FIG, panel_label

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "bench", "results")
DATA = os.path.join(HERE, "..", "figures", "data")

PANELS = [  # (experiment, family, title, x label)
    ("rings_uniform", "k8_s50_o10", r"uniform rings, $k=8$, $o_n=0.10$", "rewired fraction $p$"),
    ("rings_uniform", "k16_s100_o05", r"uniform rings, $k=16$, $o_n=0.05$", "rewired fraction $p$"),
    ("rings_hetero", "two_worlds", "two worlds", "rewired fraction $p$"),
    ("rings_hetero", "gradient10", "gradient of ten rings", "rewired fraction $p$"),
    ("rings_hetero", "sizes_only_k8", r"equal $k$, unequal sizes", "rewired fraction $p$"),
    ("lfr", "on10", r"LFR, $o_n=0.1$", r"mixing $\mu$"),
]
METHODS = [  # fixed colour per method, never cycled
    ("hc", "HyperCommon", "#2a78d6", "-"),
    ("hc_tri", "triangles only (filtered CPM)", "#2a78d6", "--"),
    ("cpm", "CPM", "#1baf7a", "-"),
    ("leiden", "Leiden", "#eb6834", "-"),
    ("angel", "ANGEL", "#6250d6", "-"),
    ("infomap", "Infomap", "#e87ba4", "-"),
    ("slpa", "SLPA", "#eda100", "-"),
]

summ = {}
for exp in {p[0] for p in PANELS}:
    s = json.load(open(os.path.join(RES, f"{exp}_summary.json")))
    for x in s["summary"]:
        summ[(exp, x["family"], x["method"], x["protocol"])] = x

fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.4), sharey=True)
rows = []
for ax, (exp, fam, title, xl), letter in zip(axes.ravel(), PANELS, "abcdef"):
    for m, lab, col, ls in METHODS:
        x = summ.get((exp, fam, m, "oracle"))
        if x is None:
            continue
        ax.fill_between(x["p"], x["lo"], x["hi"], color=col, alpha=0.15, lw=0)
        ax.plot(x["p"], x["omega"], ls, color=col, lw=1.5 if m != "hc" else 2.0, label=lab)
        rows += [dict(experiment=exp, family=fam, method=m, protocol="oracle", p=p, omega=o, lo=l, hi=h)
                 for p, o, l, h in zip(x["p"], x["omega"], x["lo"], x["hi"])]
    ax.set_title(title, fontsize=7.8)
    ax.set_xlabel(xl, fontsize=7)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    panel_label(ax, letter)
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\Omega$ (best parameter)")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=7, frameon=False, fontsize=6.6, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=1.0)
fig.savefig(os.path.join(FIG, "fig_benchmarks.pdf"))
fig.savefig(os.path.join(FIG, "fig_benchmarks.png"), dpi=200)
pd.DataFrame(rows).to_csv(os.path.join(DATA, "benchmarks_curves.csv"), index=False)
