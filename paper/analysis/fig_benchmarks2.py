"""
fig_benchmarks.pdf (label-free protocol) and fig_benchmarks_oracle.pdf (oracle protocol):
Omega against noise for representative families of the second benchmark round, mean
over the ten test graphs per noise level with 95% bootstrap bands. Also writes
figures/data/benchmarks_curves2.csv.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_real_networks import GRID, FIG, panel_label

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.join(HERE, "..", "bench", "results2")
DATA = os.path.join(HERE, "..", "figures", "data")

PANELS = [  # (experiment, family, title, x label)
    ("rings_uniform", "k8_s50_o10", r"uniform rings, $k=8$, $o_n=0.10$", "rewired fraction $p$"),
    ("rings_uniform", "k16_s100_o05", r"uniform rings, $k=16$, $o_n=0.05$", "rewired fraction $p$"),
    ("rings_hetero", "two_worlds", "two worlds", "rewired fraction $p$"),
    ("rings_hetero", "three_big_many_small", "three big, many small", "rewired fraction $p$"),
    ("rings_hetero", "gradient10", "gradient of ten rings", "rewired fraction $p$"),
    ("lfr", "on10", r"LFR, $o_n=0.1$", r"mixing $\mu$"),
]
METHODS = [  # fixed colour per method (Okabe-Ito, colour-blind safe), labels as in the tables
    ("hc", "HC", "#0072B2", "-", 2.2),
    ("hc_tri", "tri.", "#0072B2", "--", 1.2),
    ("linkcomm", "Link", "#7a7f87", "-", 1.2),
    ("cpm", "CPM", "#009E73", "-", 1.2),
    ("leiden", "Leiden", "#D55E00", "-", 1.2),
    ("leiden_cpm", "L-CPM", "#CC79A7", "-", 1.4),
    ("infomap", "Infomap", "#56B4E9", "-", 1.2),
    ("slpa", "SLPA", "#E69F00", "-", 1.2),
]
rng = np.random.default_rng(7)


def curves(protocol):
    out = []
    for exp in sorted({p[0] for p in PANELS}):
        fams = [p[1] for p in PANELS if p[0] == exp]
        d = pd.read_csv(os.path.join(R2, f"{exp}.csv"), usecols=["family", "p", "seed", "split", "method", "omega", "eq"])
        d = d[(d.split == "test") & d.family.isin(fams) & d.method.isin([m[0] for m in METHODS]) & d.omega.notna()]
        col = "omega" if protocol == "oracle" else "eq"
        sel = d.loc[d.groupby(["family", "method", "p", "seed"])[col].idxmax()]
        for (fam, m, p), g in sel.groupby(["family", "method", "p"]):
            x = g.omega.values
            b = x[rng.integers(0, len(x), (2000, len(x)))].mean(1)
            out.append(dict(experiment=exp, family=fam, method=m, protocol=protocol, p=p, omega=x.mean(),
                            lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5)))
    return pd.DataFrame(out)


def figure(C, protocol, name):
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.5), sharey=True)
    for ax, (exp, fam, title, xl), letter in zip(axes.ravel(), PANELS, "abcdef"):
        for m, lab, col, ls, lw in METHODS:
            g = C[(C.experiment == exp) & (C.family == fam) & (C.method == m)].sort_values("p")
            if g.empty:
                continue
            ax.fill_between(g.p, g.lo, g.hi, color=col, alpha=0.13, lw=0)
            ax.plot(g.p, g.omega, ls, color=col, lw=lw, label=lab)
        ax.set_title(title, fontsize=7.8)
        ax.set_xlabel(xl, fontsize=7)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        panel_label(ax, letter)
    ylab = r"$\Omega$ (label-free parameter)" if protocol == "label_free" else r"$\Omega$ (best parameter)"
    for ax in axes[:, 0]:
        ax.set_ylabel(ylab)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=8, frameon=False, fontsize=6.2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=1.0)
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    C = pd.concat([curves("label_free"), curves("oracle")], ignore_index=True)
    C.to_csv(os.path.join(DATA, "benchmarks_curves2.csv"), index=False)
    figure(C[C.protocol == "label_free"], "label_free", "fig_benchmarks")
    figure(C[C.protocol == "oracle"], "oracle", "fig_benchmarks_oracle")
    print("ok")
