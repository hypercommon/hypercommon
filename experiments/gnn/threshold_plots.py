"""
threshold_plots — omega-vs-threshold figures from the GNN omega(p, t) dataset.

Every trajectory in dataset_v1 holds a 101x101 grid: omega for each
(p, t) pair, p and t both stepping 0.01 from 0 to 1. Five runs share each
(n, z, ring, overlap) config, so every curve here is a mean over those runs
with a +-1 std band.

Putting p on the same axes rather than in separate files is what keeps the
output readable: one figure per (z, ring, overlap) with a line per p value is
198 figures, where one figure per (z, ring, overlap, run, p) would be 99,990.

Figure families (each writes into its own subdirectory):

  by_p/       one figure per (z, ring, overlap); a line per p value
              -> how the threshold response degrades as the graph is rewired
  by_ring/    one figure per (z, overlap, p); a line per ring size
              -> the ring-size comparison, holding everything else fixed
  by_z/       one figure per (ring, overlap, p); a line per z
  by_overlap/ one figure per (z, ring, p); a line per overlap %
  heatmap/    one figure per (z, ring, overlap); the full omega(p, t) surface
  summary/    t* and plateau width against p, aggregated over configs

Run from venv:
  ./.venv/Scripts/python.exe -m experiments.gnn.threshold_plots
  ./.venv/Scripts/python.exe -m experiments.gnn.threshold_plots --families by_ring,heatmap
"""

from __future__ import annotations

import argparse
import os
from functools import lru_cache

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

DATASET = os.path.join("results", "gnn", "rewire_random", "dataset_v1", "trajectories")
OUT_ROOT = os.path.join("results", "threshold_plots")

N = 400
ZS = [6, 8, 10, 12, 14, 16]
RINGS = [50, 100, 200]
OVERLAPS = list(range(0, 11))
RUNS = [1, 2, 3, 4, 5]

# p values drawn as separate lines in the by_p family. The full 101 would be
# unreadable; these span the useful range of the transition.
P_LINES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

# p values that get their own figure in the by_ring / by_z / by_overlap families.
P_SLICES = [0.0, 0.2, 0.4, 0.6]

# omega below this is treated as "no signal" when trimming the t axis.
SIGNAL_FLOOR = 0.02

# fraction of the peak that still counts as "as good as optimal"
PLATEAU_FRAC = 0.95

DPI = 130


# ---------------------------------------------------------------------------
# Data access
# ---------------------------------------------------------------------------

def trajectory_dir(z: int, ring: int, overlap_pct: int, run: int) -> str:
    return os.path.join(DATASET, f"n{N}_z{z}_r{ring}_o{overlap_pct:02d}_run{run}")


@lru_cache(maxsize=None)
def load_config(z: int, ring: int, overlap_pct: int):
    """Mean and std of omega(p, t) over the available runs.

    Returns (mean, std, p_grid, t_grid, n_runs) or None when nothing is on disk.
    """
    grids = []
    p_grid = t_grid = None

    for run in RUNS:
        path = os.path.join(trajectory_dir(z, ring, overlap_pct, run), "omega_grid.npz")
        if not os.path.exists(path):
            continue
        data = np.load(path)
        grids.append(data["omega"])
        if p_grid is None:
            p_grid, t_grid = data["p_grid"], data["t_grid"]

    if not grids:
        return None

    stack = np.stack(grids)
    return stack.mean(axis=0), stack.std(axis=0), p_grid, t_grid, len(grids)


def t_limit(*mean_grids) -> float:
    """Largest t worth drawing — beyond the signal, omega is pinned near zero.

    At t=0 every pair passes the predicate and the graph collapses into one
    community; at t=1 nothing passes and there are no communities. Both ends
    score ~0, so plotting the full [0, 1] leaves most of the axis empty.
    """
    highest = 0.0
    for mean in mean_grids:
        if mean is None:
            continue
        columns = np.where(mean.max(axis=0) > SIGNAL_FLOOR)[0]
        if len(columns):
            highest = max(highest, columns.max() / (mean.shape[1] - 1))
    return min(1.0, max(0.25, highest + 0.06))


def plateau_width(row: np.ndarray, t_grid: np.ndarray) -> float:
    """Width in t of the region scoring at least PLATEAU_FRAC of the peak."""
    peak = np.nanmax(row)
    if not np.isfinite(peak) or peak <= SIGNAL_FLOOR:
        return 0.0
    inside = np.where(row >= PLATEAU_FRAC * peak)[0]
    if len(inside) == 0:
        return 0.0
    step = t_grid[1] - t_grid[0]
    return float((inside.max() - inside.min() + 1) * step)


def save(fig, *parts: str) -> str:
    path = os.path.join(OUT_ROOT, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def style_axes(ax, t_max: float, title: str) -> None:
    ax.set_xlabel("threshold  t")
    ax.set_ylabel(r"$\omega$  (vs ground truth)")
    ax.set_title(title)
    ax.set_xlim(0, t_max)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.3, linewidth=0.6)
    ax.axhline(0, color="black", linewidth=0.8, alpha=0.5)


# ---------------------------------------------------------------------------
# Families
# ---------------------------------------------------------------------------

def family_by_p() -> int:
    """One figure per (z, ring, overlap); a line per p value."""
    made = 0
    for z in ZS:
        for ring in RINGS:
            for overlap in OVERLAPS:
                loaded = load_config(z, ring, overlap)
                if loaded is None:
                    continue
                mean, std, p_grid, t_grid, n_runs = loaded
                t_max = t_limit(mean)

                fig, ax = plt.subplots(figsize=(7.2, 4.6))
                colors = plt.cm.viridis(np.linspace(0, 0.92, len(P_LINES)))

                for color, p in zip(colors, P_LINES):
                    i = int(round(p * (len(p_grid) - 1)))
                    ax.plot(t_grid, mean[i], color=color, linewidth=1.6, label=f"p={p:.1f}")
                    ax.fill_between(
                        t_grid, mean[i] - std[i], mean[i] + std[i],
                        color=color, alpha=0.15, linewidth=0,
                    )

                style_axes(
                    ax, t_max,
                    f"n={N}, z={z}, ring={ring}, overlap={overlap}%"
                    f"   (mean of {n_runs} runs $\\pm$1 std)",
                )
                ax.legend(fontsize=8, ncol=2, framealpha=0.9)
                save(fig, "by_p", f"z{z}_ring{ring}", f"overlap{overlap:02d}.png")
                made += 1
    return made


def family_by_ring() -> int:
    """One figure per (z, overlap, p); a line per ring size."""
    made = 0
    for z in ZS:
        for overlap in OVERLAPS:
            loaded = {r: load_config(z, r, overlap) for r in RINGS}
            if not any(v is not None for v in loaded.values()):
                continue
            t_max = t_limit(*[v[0] if v else None for v in loaded.values()])

            for p in P_SLICES:
                fig, ax = plt.subplots(figsize=(7.2, 4.6))
                colors = plt.cm.plasma(np.linspace(0.1, 0.8, len(RINGS)))

                for color, ring in zip(colors, RINGS):
                    if loaded[ring] is None:
                        continue
                    mean, std, p_grid, t_grid, n_runs = loaded[ring]
                    i = int(round(p * (len(p_grid) - 1)))
                    rings_count = N // ring
                    ax.plot(
                        t_grid, mean[i], color=color, linewidth=1.8,
                        label=f"ring={ring}  ({rings_count} rings)",
                    )
                    ax.fill_between(
                        t_grid, mean[i] - std[i], mean[i] + std[i],
                        color=color, alpha=0.15, linewidth=0,
                    )

                style_axes(ax, t_max, f"n={N}, z={z}, overlap={overlap}%, p={p:.2f}")
                ax.legend(fontsize=9, framealpha=0.9)
                save(fig, "by_ring", f"z{z}_overlap{overlap:02d}", f"p{p:.2f}.png")
                made += 1
    return made


def family_by_z() -> int:
    """One figure per (ring, overlap, p); a line per z."""
    made = 0
    for ring in RINGS:
        for overlap in OVERLAPS:
            loaded = {z: load_config(z, ring, overlap) for z in ZS}
            if not any(v is not None for v in loaded.values()):
                continue
            t_max = t_limit(*[v[0] if v else None for v in loaded.values()])

            for p in P_SLICES:
                fig, ax = plt.subplots(figsize=(7.2, 4.6))
                colors = plt.cm.viridis(np.linspace(0, 0.92, len(ZS)))

                for color, z in zip(colors, ZS):
                    if loaded[z] is None:
                        continue
                    mean, std, p_grid, t_grid, n_runs = loaded[z]
                    i = int(round(p * (len(p_grid) - 1)))
                    ax.plot(t_grid, mean[i], color=color, linewidth=1.7, label=f"z={z}")
                    ax.fill_between(
                        t_grid, mean[i] - std[i], mean[i] + std[i],
                        color=color, alpha=0.13, linewidth=0,
                    )

                style_axes(ax, t_max, f"n={N}, ring={ring}, overlap={overlap}%, p={p:.2f}")
                ax.legend(fontsize=9, ncol=2, framealpha=0.9)
                save(fig, "by_z", f"ring{ring}_overlap{overlap:02d}", f"p{p:.2f}.png")
                made += 1
    return made


def family_by_overlap() -> int:
    """One figure per (z, ring, p); a line per overlap %."""
    made = 0
    for z in ZS:
        for ring in RINGS:
            loaded = {o: load_config(z, ring, o) for o in OVERLAPS}
            if not any(v is not None for v in loaded.values()):
                continue
            t_max = t_limit(*[v[0] if v else None for v in loaded.values()])

            for p in P_SLICES:
                fig, ax = plt.subplots(figsize=(7.2, 4.6))
                colors = plt.cm.coolwarm(np.linspace(0, 1, len(OVERLAPS)))

                for color, overlap in zip(colors, OVERLAPS):
                    if loaded[overlap] is None:
                        continue
                    mean, std, p_grid, t_grid, n_runs = loaded[overlap]
                    i = int(round(p * (len(p_grid) - 1)))
                    ax.plot(t_grid, mean[i], color=color, linewidth=1.5, label=f"{overlap}%")

                style_axes(ax, t_max, f"n={N}, z={z}, ring={ring}, p={p:.2f}")
                ax.legend(fontsize=8, ncol=3, title="overlap", framealpha=0.9)
                save(fig, "by_overlap", f"z{z}_ring{ring}", f"p{p:.2f}.png")
                made += 1
    return made


def family_heatmap() -> int:
    """One figure per (z, ring, overlap): the whole omega(p, t) surface."""
    made = 0
    for z in ZS:
        for ring in RINGS:
            for overlap in OVERLAPS:
                loaded = load_config(z, ring, overlap)
                if loaded is None:
                    continue
                mean, std, p_grid, t_grid, n_runs = loaded
                t_max = t_limit(mean)
                keep = t_grid <= t_max

                fig, ax = plt.subplots(figsize=(7.0, 4.8))
                im = ax.imshow(
                    mean[:, keep],
                    aspect="auto", origin="lower", cmap="magma",
                    vmin=0, vmax=1,
                    extent=[t_grid[keep].min(), t_grid[keep].max(), p_grid.min(), p_grid.max()],
                )
                # t* per p, drawn over the surface
                best = np.array([t_grid[int(np.nanargmax(row))] for row in mean])
                ax.plot(best, p_grid, color="cyan", linewidth=1.2, alpha=0.9, label="$t^*(p)$")

                ax.set_xlabel("threshold  t")
                ax.set_ylabel("rewiring probability  p")
                ax.set_title(
                    f"n={N}, z={z}, ring={ring}, overlap={overlap}%"
                    f"   (mean of {n_runs} runs)"
                )
                ax.legend(fontsize=8, loc="upper right")
                fig.colorbar(im, ax=ax, label=r"$\omega$")
                save(fig, "heatmap", f"z{z}_ring{ring}", f"overlap{overlap:02d}.png")
                made += 1
    return made


def family_summary() -> int:
    """t* and plateau width against p, one panel per ring size."""
    made = 0

    # Ring size is encoded as line style so a single legend can carry both
    # dimensions without one entry per (z, ring) combination.
    ring_styles = dict(zip(RINGS, ["-", "--", ":"]))

    for overlap in [0, 5, 10]:
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
        colors = dict(zip(ZS, plt.cm.viridis(np.linspace(0, 0.92, len(ZS)))))

        for ring in RINGS:
            for z in ZS:
                loaded = load_config(z, ring, overlap)
                if loaded is None:
                    continue
                mean, std, p_grid, t_grid, n_runs = loaded

                best, width = [], []
                for row in mean:
                    if np.nanmax(row) <= SIGNAL_FLOOR:
                        best.append(np.nan)
                        width.append(np.nan)
                    else:
                        best.append(t_grid[int(np.nanargmax(row))])
                        width.append(plateau_width(row, t_grid))

                for ax, series in ((axes[0], best), (axes[1], width)):
                    ax.plot(
                        p_grid, series,
                        color=colors[z], linestyle=ring_styles[ring],
                        linewidth=1.2, alpha=0.85,
                    )

        # one legend entry per z (colour) and one per ring size (line style)
        handles = [
            plt.Line2D([], [], color=colors[z], linewidth=1.6, label=f"z={z}") for z in ZS
        ] + [
            plt.Line2D([], [], color="grey", linestyle=ring_styles[r], linewidth=1.6,
                       label=f"ring={r}")
            for r in RINGS
        ]

        axes[0].set_xlabel("rewiring probability  p")
        axes[0].set_ylabel("$t^*$  (argmax $\\omega$)")
        axes[0].set_title(f"best threshold vs p   (overlap={overlap}%, all ring sizes)")
        axes[1].set_xlabel("rewiring probability  p")
        axes[1].set_ylabel(f"width of $\\omega \\geq$ {PLATEAU_FRAC:.2f}$\\cdot$peak")
        axes[1].set_title(f"plateau width vs p   (overlap={overlap}%)")
        for ax in axes:
            ax.grid(alpha=0.3, linewidth=0.6)
            ax.legend(handles=handles, fontsize=7, ncol=3, framealpha=0.9)

        save(fig, "summary", f"tstar_and_plateau_overlap{overlap:02d}.png")
        made += 1

    return made


FAMILIES = {
    "by_p": family_by_p,
    "by_ring": family_by_ring,
    "by_z": family_by_z,
    "by_overlap": family_by_overlap,
    "heatmap": family_heatmap,
    "summary": family_summary,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--families",
        default=",".join(FAMILIES),
        help="comma-separated subset of: " + ", ".join(FAMILIES),
    )
    args = parser.parse_args()

    chosen = [f.strip() for f in args.families.split(",") if f.strip()]
    unknown = [f for f in chosen if f not in FAMILIES]
    if unknown:
        raise SystemExit(f"unknown families: {unknown}. available: {list(FAMILIES)}")

    os.makedirs(OUT_ROOT, exist_ok=True)
    total = 0
    for name in chosen:
        count = FAMILIES[name]()
        total += count
        print(f"  {name:12} {count:5} figures")

    print(f"\n{total} figures written under {OUT_ROOT}")


if __name__ == "__main__":
    main()
