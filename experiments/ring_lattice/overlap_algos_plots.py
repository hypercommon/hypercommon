"""Plot omega vs p for all algos, one figure per overlap level."""

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt


ALGOS = [
    "hypercommon",
    "leiden",
    "label_propagation",
    "walktrap",
    "slpa",
    "demon",
    "angel",
]


def overlap_label(overlap: float) -> str:
    pct = round(overlap * 100)
    return f"o{pct}"


def infer_config(run_dir: str) -> str:
    """Best-effort config string for plot titles. Reads parent dir name + first
    algo's progress.csv. Falls back to a generic label."""
    parent = os.path.basename(os.path.dirname(run_dir))
    # parse e.g. 'overlap_algos_z8_r50' -> 'z=8, ring=50'
    parts = parent.replace("overlap_algos", "").strip("_").split("_")
    pieces = []
    for p in parts:
        if p.startswith("z"):
            pieces.append(f"z={p[1:]}")
        elif p.startswith("r"):
            pieces.append(f"ring={p[1:]}")
        elif p.startswith("n"):
            pieces.append(f"n={p[1:]}")
    if pieces:
        return ", ".join(pieces)
    return parent


def plot_run(run_dir: str):
    """Read each algo's avg.xlsx and produce one plot per overlap."""
    plots_dir = os.path.join(run_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    config_label = infer_config(run_dir)

    # Load every algo's avg table
    by_algo: dict[str, pd.DataFrame] = {}
    for algo in ALGOS:
        path = os.path.join(run_dir, algo, "avg.xlsx")
        if not os.path.exists(path):
            print(f"missing: {path}")
            continue
        by_algo[algo] = pd.read_excel(path, sheet_name="avg")

    if not by_algo:
        raise FileNotFoundError(f"no avg.xlsx files found under {run_dir}")

    overlaps = sorted(set().union(*[set(df["overlap"].unique()) for df in by_algo.values()]))
    print(f"overlaps found: {overlaps}")
    print(f"algos found: {list(by_algo.keys())}")

    # tab10 has 10 distinct colours — plenty for 7 algos
    colours = plt.cm.tab10.colors
    algo_colour = {algo: colours[i] for i, algo in enumerate(ALGOS)}

    for ov in overlaps:
        plt.figure(figsize=(11, 6))
        for algo, df in by_algo.items():
            sub = df[df["overlap"] == ov].sort_values("p")
            if sub.empty:
                continue
            plt.plot(sub["p"], sub["omega_avg"], color=algo_colour[algo], linewidth=1.6, label=algo)

        plt.xlabel("p (fraction rewired)")
        plt.ylabel("omega (avg over 10 runs)")
        plt.title(f"omega vs p — all algos — overlap = {round(ov*100)}%  ({config_label})")
        plt.xlim(0, 1)
        plt.ylim(-0.05, 1.05)
        plt.xticks([i / 10 for i in range(11)])
        plt.yticks([i / 10 for i in range(11)])
        plt.grid(True, alpha=0.25)
        plt.legend(loc="best", frameon=True)
        plt.tight_layout()

        path = os.path.join(plots_dir, f"omega_vs_p_{overlap_label(ov)}.png")
        plt.savefig(path, dpi=200)
        plt.close()
        print(f"saved {path}")


def latest_run_dir(out_dir: str = "results/overlap_algos") -> str:
    candidates = sorted(glob.glob(os.path.join(out_dir, "run_*")))
    if not candidates:
        raise FileNotFoundError(f"no runs found under {out_dir}")
    return candidates[-1]


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else latest_run_dir()
    print(f"target: {target}")
    plot_run(target)
