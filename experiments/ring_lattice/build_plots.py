"""
Generate all summary plots: for each combo of 3 fixed params, vary the 4th on x-axis.
3 plot types per combo: t_best only, p_crit only, both on same axes.
Output: results/pcrit_experiment/final_result/plots/x=<var>/
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import product

SUMMARY = "results/pcrit_experiment/final_result/summary.xlsx"
PLOTS_ROOT = "results/pcrit_experiment/final_result/plots"

COLOR_T    = "#2196F3"
COLOR_P    = "#E91E63"
MARKER     = "o"


def save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _fmt_tick(v):
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _set_data_xticks(ax, xs):
    ax.set_xticks(list(xs))
    ax.set_xticklabels([_fmt_tick(v) for v in xs])


def _set_data_yticks(ax, *value_lists):
    ticks = sorted({round(v, 2) for vs in value_lists for v in vs if v is not None})
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:.2f}" for t in ticks])
    if ticks:
        pad = max(0.005, (ticks[-1] - ticks[0]) * 0.05)
        ax.set_ylim(ticks[0] - pad, ticks[-1] + pad)


def _fig_with_table(rows_data, row_labels, col_labels):
    fig, (ax_plot, ax_tbl) = plt.subplots(
        2, 1, figsize=(10, 7),
        gridspec_kw={"height_ratios": [4, 1]},
    )
    col_text = [_fmt_tick(c) for c in col_labels]
    cell_text = [[f"{v:.2f}" if v is not None else "N/A" for v in row] for row in rows_data]
    ax_tbl.axis("off")
    tbl = ax_tbl.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=col_text,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.auto_set_column_width(list(range(len(col_text))))
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.5)
        if r == 0 or c == -1:
            cell.set_facecolor("#d0d8e8")
    return fig, ax_plot


def plot_single(xs, ys, xlabel, ylabel, color, title, path):
    fig, ax = _fig_with_table([ys], [ylabel], xs)
    ax.plot(xs, ys, marker=MARKER, color=color, linewidth=1.8, markersize=7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    _set_data_xticks(ax, xs)
    _set_data_yticks(ax, ys)
    ax.grid(True, alpha=0.3, linestyle="--")
    save(fig, path)


def plot_both(xs, yt, yp, xlabel, title, path):
    fig, ax = _fig_with_table([yp, yt], ["p_crit", "t_best"], xs)
    ax.plot(xs, yp, marker="o", color="steelblue", linewidth=1.8, markersize=7, label="p_crit")
    ax.plot(xs, yt, marker="s", linestyle="--", color="darkorange",
            linewidth=1.8, markersize=7, label="t_best")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("value")
    ax.set_title(title, fontsize=11)
    _set_data_xticks(ax, xs)
    _set_data_yticks(ax, yt, yp)
    ax.legend(loc="best", fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle="--")
    save(fig, path)


def make_plots_for_xvar(df, xvar, fixed_vars):
    folder = os.path.join(PLOTS_ROOT, f"x={xvar}")
    xvals_all = sorted(df[xvar].unique())

    fixed_combos = list(product(*[sorted(df[v].unique()) for v in fixed_vars]))
    count = 0

    for combo in fixed_combos:
        mask = pd.Series([True] * len(df), index=df.index)
        for v, val in zip(fixed_vars, combo):
            mask &= df[v] == val
        sub = df[mask].sort_values(xvar)

        if len(sub) < 2:
            continue

        xs = sub[xvar].tolist()
        yt = sub["t_best"].tolist()
        yp = sub["p_crit"].tolist()

        fixed_str = "_".join(
            f"{v}{val:.2f}".rstrip("0").rstrip(".") if isinstance(val, float) else f"{v}{val}"
            for v, val in zip(fixed_vars, combo)
        )
        title_fixed = "  ".join(
            f"{v}={val:.2f}".rstrip("0").rstrip(".") if isinstance(val, float) else f"{v}={val}"
            for v, val in zip(fixed_vars, combo)
        )
        xlabel = xvar

        base = os.path.join(folder, fixed_str)

        plot_single(xs, yt, xlabel, "t_best", COLOR_T,
                    f"t_best  |  {title_fixed}", f"{base}_tbest.png")
        plot_single(xs, yp, xlabel, "p_crit", COLOR_P,
                    f"p_crit  |  {title_fixed}", f"{base}_pcrit.png")
        plot_both(xs, yt, yp, xlabel,
                  f"t_best & p_crit  |  {title_fixed}", f"{base}_both.png")
        count += 1

    return count * 3


def main():
    df = pd.read_excel(SUMMARY)

    axes = [
        ("overlap",   ["n", "z", "ring_size"]),
        ("n",         ["z", "ring_size", "overlap"]),
        ("z",         ["n", "ring_size", "overlap"]),
        ("ring_size", ["n", "z", "overlap"]),
    ]

    total = 0
    for xvar, fixed_vars in axes:
        n_plots = make_plots_for_xvar(df, xvar, fixed_vars)
        print(f"x={xvar:<12}  {n_plots} plots")
        total += n_plots

    print(f"\nTotal: {total} plots -> {PLOTS_ROOT}")


if __name__ == "__main__":
    main()
