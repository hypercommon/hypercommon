"""
shapes_omega_vs_p - omega against rewiring probability, from a sweep's records.csv.

Reads only what a sweep already wrote; runs no algorithms, so figures regenerate
in seconds and a partially complete sweep plots fine.

Two families, both keyed on (shape, overlap):

  per_algo/<spec>/<algo>/overlap<pct>.png    one algorithm per figure
  combined/<spec>/overlap<pct>.png           all algorithms on one axes

Each curve is the mean over the sweep's runs. The shape is captioned by the
(ring size, degree) pairs it is built from, e.g. "(200,8)x4" for four rings of
200 nodes at degree 8.

Algorithm colours match the original n=2000 comparison figures under
results/overlap_algos/: those took tab10 indexed by position in the seven-algo
list, so hypercommon is blue and leiden orange. The hex values are pinned below
rather than recomputed, because the list has since dropped slpa and
label_propagation and positional indexing would silently recolour everything.

Run from venv:
  ./.venv/Scripts/python.exe -m plots.shapes_omega_vs_p
  ./.venv/Scripts/python.exe -m plots.shapes_omega_vs_p --sweep shapes_combined
  ./.venv/Scripts/python.exe -m plots.shapes_omega_vs_p --run-dir results/shapes_sizes/run_20260729_174837
  ./.venv/Scripts/python.exe -m plots.shapes_omega_vs_p --family combined
"""

from __future__ import annotations

import argparse
import glob
import os
from collections import Counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

# Frozen so these figures stay comparable with the earlier n=2000 ones.
ALGO_COLORS = {
    "hypercommon": "#1f77b4",
    "leiden": "#ff7f0e",
    "label_propagation": "#2ca02c",
    "walktrap": "#d62728",
    "slpa": "#9467bd",
    "demon": "#8c564b",
    "angel": "#e377c2",
}

# Legend and drawing order; anything unlisted is appended alphabetically.
ALGO_ORDER = [
    "hypercommon", "leiden", "label_propagation",
    "walktrap", "slpa", "demon", "angel",
]

SWEEPS = ["shapes_algos", "shapes_sizes", "shapes_combined"]

DPI = 200


def latest_run_dir(sweep: str) -> str:
    runs = sorted(glob.glob(os.path.join("results", sweep, "run_*")))
    if not runs:
        raise SystemExit(f"no runs found under results/{sweep}")
    return runs[-1]


def load(run_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Concatenate the part files, resolving reruns.

    A sweep that was paused, resumed or had shapes launched by hand writes
    several parts per shape (records_<shape>.csv, _r2, _x). A unit that was
    interrupted leaves partial rows behind that the rerun then wrote again, so
    later parts win on (run_id, p, algo).
    """
    parts = glob.glob(os.path.join(run_dir, "records_*.csv"))
    parts = [p for p in parts if os.path.getsize(p) > 200]
    if not parts:
        raise SystemExit(f"no records_*.csv in {run_dir}")

    parts.sort(key=lambda p: (("_r2" in p) or ("_x" in p), p))
    frame = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    frame = frame.drop_duplicates(subset=["run_id", "p", "algo"], keep="last")
    frame = frame.dropna(subset=["omega"])

    shapes_path = os.path.join(run_dir, "shapes.csv")
    shapes = pd.read_csv(shapes_path) if os.path.exists(shapes_path) else pd.DataFrame()

    return frame, shapes


def averaged(records: pd.DataFrame) -> pd.DataFrame:
    """Mean over runs for each (shape, overlap, algo, p)."""
    grouped = records.groupby(["shape", "overlap_pct", "algo", "p"], as_index=False)
    return grouped.agg(omega_mean=("omega", "mean"))


def ordered_algos(present) -> list[str]:
    present = set(present)
    known = [a for a in ALGO_ORDER if a in present]
    return known + sorted(present - set(known))


def shape_terms(shapes: pd.DataFrame, shape: str) -> list[tuple[int, int, int]]:
    """The (ring size, z, count) terms a shape is built from, largest first.

    Rings sharing a (size, z) collapse into one term, so (200, 8, 4) is four
    rings of 200 nodes at degree 8.
    """
    if shapes.empty or shape not in set(shapes["shape"]):
        return []

    row = shapes[shapes["shape"] == shape].iloc[0]
    sizes = [int(x) for x in str(row["sizes"]).split()]
    zs = [int(x) for x in str(row["zs"]).split()]

    counts = Counter(zip(sizes, zs))
    return [(size, z, count)
            for (size, z), count in sorted(counts.items(), key=lambda kv: (-kv[0][0], -kv[0][1]))]


def shape_spec(shapes: pd.DataFrame, shape: str) -> str:
    """The shape spec as it appears in a figure title.

    The shape's own name is deliberately unused — it is an arbitrary label, and
    this says the same thing exactly.
    """
    return ", ".join(f"({size},{z})*{count}"
                     for size, z, count in shape_terms(shapes, shape))


def shape_dirname(shapes: pd.DataFrame, shape: str) -> str:
    """The same spec as a directory name.

    Windows forbids '*' in a path and dislikes spaces, so the count separator
    becomes 'x' and terms join on '_'. Falls back to the shape's name if the
    spec is unavailable, so a directory is always produced.
    """
    terms = shape_terms(shapes, shape)
    if not terms:
        return shape
    return "_".join(f"({size},{z})x{count}" for size, z, count in terms)


def style(ax, title: str) -> None:
    ax.set_xlabel("p  (fraction of edges rewired)")
    ax.set_ylabel(r"$\omega$")
    # One centred line: what the run is, then what the graph is made of. Long
    # specs ease the font down rather than wrap, so the header never grows.
    size = 13.0 if len(title) <= 60 else (11.0 if len(title) <= 100 else 9.5)
    ax.set_title(title, fontsize=size, pad=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.03)
    ax.set_xticks([i / 10 for i in range(11)])
    ax.set_yticks([i / 10 for i in range(11)])
    ax.grid(True, alpha=0.18, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, *parts: str) -> None:
    path = os.path.join(*parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def family_per_algo(avg: pd.DataFrame, shapes: pd.DataFrame, plots_dir: str) -> int:
    """One figure per (shape, overlap, algo)."""
    made = 0
    for (shape, overlap, algo), frame in avg.groupby(["shape", "overlap_pct", "algo"]):
        frame = frame.sort_values("p")

        fig, ax = plt.subplots(figsize=(9, 5.2))
        ax.plot(frame["p"], frame["omega_mean"],
                color=ALGO_COLORS.get(algo, "#333333"), linewidth=1.5, label=algo)

        style(ax, f"{shape_spec(shapes, shape)}   |   {algo}   |   overlap {overlap}%")
        ax.legend(loc="best", frameon=False)
        save(fig, plots_dir, "per_algo", shape_dirname(shapes, shape), algo,
             f"overlap{overlap:02d}.png")
        made += 1
    return made


def family_combined(avg: pd.DataFrame, shapes: pd.DataFrame, plots_dir: str) -> int:
    """One figure per (shape, overlap), every algorithm together."""
    made = 0
    for (shape, overlap), frame in avg.groupby(["shape", "overlap_pct"]):
        fig, ax = plt.subplots(figsize=(10.5, 5.8))

        for algo in ordered_algos(frame["algo"].unique()):
            series = frame[frame["algo"] == algo].sort_values("p")
            # Every algorithm gets the same weight: emphasising the subject of
            # the comparison would put a thumb on the scale. hypercommon only
            # draws last so it stays visible where curves cross.
            ax.plot(series["p"], series["omega_mean"],
                    color=ALGO_COLORS.get(algo, "#333333"),
                    linewidth=1.5,
                    zorder=3 if algo == "hypercommon" else 2,
                    label=algo)

        style(ax, f"{shape_spec(shapes, shape)}   |   overlap {overlap}%")
        ax.legend(loc="best", frameon=False)
        save(fig, plots_dir, "combined", shape_dirname(shapes, shape),
             f"overlap{overlap:02d}.png")
        made += 1
    return made


FAMILIES = {"per_algo": family_per_algo, "combined": family_combined}


def plot_run(run_dir: str, families: list[str]) -> None:
    records, shapes = load(run_dir)
    avg = averaged(records)
    plots_dir = os.path.join(run_dir, "plots")

    print(f"run: {run_dir}")
    print(f"  {len(records):,} rows | {records['shape'].nunique()} shapes | "
          f"{records['algo'].nunique()} algos | overlaps {sorted(records.overlap_pct.unique())}")

    for name in families:
        count = FAMILIES[name](avg, shapes, plots_dir)
        print(f"  {name:10} {count:4} figures")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=None,
                        help="a single run directory; overrides --sweep")
    parser.add_argument("--sweep", default=None, choices=SWEEPS,
                        help="plot the latest run of this sweep (default: all three)")
    parser.add_argument("--family", default=",".join(FAMILIES),
                        help="comma-separated subset of: " + ", ".join(FAMILIES))
    args = parser.parse_args()

    families = [f.strip() for f in args.family.split(",") if f.strip()]
    unknown = [f for f in families if f not in FAMILIES]
    if unknown:
        raise SystemExit(f"unknown families: {unknown}. available: {list(FAMILIES)}")

    if args.run_dir:
        run_dirs = [args.run_dir]
    elif args.sweep:
        run_dirs = [latest_run_dir(args.sweep)]
    else:
        run_dirs = [latest_run_dir(s) for s in SWEEPS
                    if glob.glob(os.path.join("results", s, "run_*"))]

    for run_dir in run_dirs:
        plot_run(run_dir, families)


if __name__ == "__main__":
    main()
