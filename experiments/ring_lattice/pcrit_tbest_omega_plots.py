"""Plot averaged t_argmax(p) and omega_at_argmax(p) over the 3 runs."""

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt


def plot_run(run_dir: str):
    run_subdirs = sorted(glob.glob(os.path.join(run_dir, "run_*")))
    if not run_subdirs:
        raise FileNotFoundError(f"no run_* subdirs in {run_dir}")

    print(f"Reading {len(run_subdirs)} run subdirs:")
    frames = []
    for d in run_subdirs:
        path = os.path.join(d, "argmax.csv")
        df = pd.read_csv(path)
        df["run"] = os.path.basename(d)
        frames.append(df)
        print(f"  {d}  ({len(df)} rows)")

    full = pd.concat(frames, ignore_index=True)
    avg = full.groupby("p", as_index=False).agg(
        t_argmax_avg=("t_argmax", "mean"),
        t_argmax_std=("t_argmax", "std"),
        omega_avg=("omega_at_argmax", "mean"),
        omega_std=("omega_at_argmax", "std"),
    )

    plots_dir = os.path.join(run_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # ---- Plot 1: t_best vs p ----
    plt.figure(figsize=(11, 6))
    plt.plot(avg["p"], avg["t_argmax_avg"], color="tab:blue", linewidth=2)
    plt.xlabel("p (fraction rewired)")
    plt.ylabel("t_argmax (best threshold)")
    plt.title("Best Jaccard threshold vs rewiring probability  (n=2000, z=16, ring=100, overlap=0)")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xticks([i / 10 for i in range(11)])
    plt.yticks([i / 10 for i in range(11)])
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    p1 = os.path.join(plots_dir, "t_argmax_vs_p.png")
    plt.savefig(p1, dpi=200)
    plt.close()
    print(f"saved {p1}")

    # ---- Plot 2: omega_at_argmax vs p ----
    plt.figure(figsize=(11, 6))
    plt.plot(avg["p"], avg["omega_avg"], color="tab:green", linewidth=2)
    plt.xlabel("p (fraction rewired)")
    plt.ylabel("omega at best t (recognition success)")
    plt.title("Best omega vs rewiring probability  (n=2000, z=16, ring=100, overlap=0)")
    plt.xlim(0, 1)
    plt.ylim(-0.05, 1.05)
    plt.xticks([i / 10 for i in range(11)])
    plt.yticks([i / 10 for i in range(11)])
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    p2 = os.path.join(plots_dir, "omega_at_argmax_vs_p.png")
    plt.savefig(p2, dpi=200)
    plt.close()
    print(f"saved {p2}")

    # ---- Plot 3: comparison vs pcrit (single fixed t_best) ----
    # final_result supersedes the per-run CSVs; it carries an extra overlap column,
    # and this comparison is against the zero-overlap sweep.
    pcrit_xlsx = "results/pcrit_experiment/final_result/p_sweep.xlsx"
    if os.path.exists(pcrit_xlsx):
        pcrit = pd.read_excel(pcrit_xlsx)
        pcrit = pcrit[(pcrit["n"] == 2000) & (pcrit["z"] == 16) & (pcrit["ring_size"] == 100)]
        if "overlap" in pcrit.columns:
            pcrit = pcrit[pcrit["overlap"] == 0.0]
        pcrit_tbest = pcrit["t_best"].iloc[0]

        plt.figure(figsize=(11, 6))
        plt.plot(avg["p"], avg["omega_avg"], color="tab:green", linewidth=2,
                 label="envelope (per-p t_argmax)")
        plt.plot(pcrit["p"], pcrit["score"], color="tab:red", linewidth=2,
                 label=f"pcrit (fixed t={pcrit_tbest})")
        plt.xlabel("p (fraction rewired)")
        plt.ylabel("omega")
        plt.title("Envelope vs pcrit single-t  (n=2000, z=16, ring=100, overlap=0)")
        plt.xlim(0, 1)
        plt.ylim(-0.05, 1.05)
        plt.xticks([i / 10 for i in range(11)])
        plt.yticks([i / 10 for i in range(11)])
        plt.grid(True, alpha=0.25)
        plt.legend()
        plt.tight_layout()
        p3 = os.path.join(plots_dir, "envelope_vs_pcrit.png")
        plt.savefig(p3, dpi=200)
        plt.close()
        print(f"saved {p3}")
    else:
        print(f"pcrit comparison file not found: {pcrit_csv}")

    # also save the averaged table for reference
    avg_path = os.path.join(plots_dir, "averaged.csv")
    avg.to_csv(avg_path, index=False)
    print(f"saved {avg_path}")


def latest_run_dir(out_dir: str = "results/pcrit_tbest_omega") -> str:
    candidates = sorted(glob.glob(os.path.join(out_dir, "run_*")))
    if not candidates:
        raise FileNotFoundError(f"no runs found under {out_dir}")
    return candidates[-1]


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else latest_run_dir()
    print(f"target: {target}")
    plot_run(target)
