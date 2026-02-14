import os
import random
import time
from fractions import Fraction
from itertools import combinations

import matplotlib.pyplot as plt

from generators.ring_lattice import ring_lattice
from utils.rewiring import rewire_step
from metrics.omega import omega_index, build_pair_counts
from predicates.jaccard import closed_neighborhood_jaccard_predicate
from hypercommon.algorithm import get_communities


def ring_ground_truth(n: int, rings: int) -> list[set[int]]:
    n_per = n // rings
    return [set(range(r * n_per, (r + 1) * n_per)) for r in range(rings)]


def safe_sheet_num(x: float) -> str:
    return f"{x:.6f}".rstrip("0").rstrip(".").replace(".", "p")


def ring_lattice_edge_count(n: int, z: int, rings: int) -> int:
    # For disjoint rings: each ring has n_per nodes, z-regular undirected -> edges = n_per*z/2
    n_per = n // rings
    return rings * (n_per * z // 2)


def validate_rewiring_plan(M0: int, p_step: float):
    p = Fraction(str(p_step)).limit_denominator()

    # steps must be integer: p_step = 1/steps
    if p.numerator != 1:
        raise ValueError(f"p_step must be of the form 1/steps (e.g. 0.01). Got {p_step} = {p}")

    steps = p.denominator  # because p_step = 1/steps

    # k_step must be integer: M0 / steps
    if M0 % steps != 0:
        raise ValueError(
            f"Illegal p_step={p_step}: M0={M0} edges is not divisible by steps={steps}. "
            f"Choose p_step so that steps divides M0."
        )

    k_step = M0 // steps
    return steps, k_step


def run_experiment(
    configs,
    avg_times: int = 10,
    p_step: float = 0.01,
    seed: int = 42,
    out_dir: str = "results/Rewiring VS True",
):
    if not (0.0 < p_step <= 1.0):
        raise ValueError("p_step must be in (0, 1]")

    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(seed)

    t_global0 = time.perf_counter()
    total_configs = len(configs)

    for ci, (n, z, rings) in enumerate(configs, start=1):
        t_cfg0 = time.perf_counter()
        print(f"\n[CONFIG {ci}/{total_configs}] n={n}, z={z}, rings={rings}  (start)")

        # ground truth
        truth = ring_ground_truth(n, rings)
        gt_pair_counts = build_pair_counts(truth)

        # omega helpers (no need for G_ref)
        total_pairs = n * (n - 1) // 2

        # lattice edges count
        M0 = ring_lattice_edge_count(n=n, z=z, rings=rings)

        steps, k_step = validate_rewiring_plan(M0, p_step)
        xs = [s / steps for s in range(steps + 1)]

        thresholds = [0.127, 0.128, 0.129, 0.13, 0.131, 0.132, 0.133]
        total_thresholds = len(thresholds)

        print(f"[CONFIG {ci}] M0={M0}, k_step={k_step}, steps={steps}, thresholds={total_thresholds}")

        for ti, threshold in enumerate(thresholds, start=1):
            t_thr0 = time.perf_counter()
            print(f"  [THR {ti}/{total_thresholds}] t={threshold:.3f} (start)")

            sum_scores = [0.0] * (steps + 1)

            for run_i in range(1, avg_times + 1):
                t_run0 = time.perf_counter()

                G = ring_lattice(n=n, z=z, rings=rings)

                edge_stack = list(G.edges())

                rng.shuffle(edge_stack)

                for s in range(steps):
                    t_step0 = time.perf_counter()

                    pred = get_communities(G, closed_neighborhood_jaccard_predicate(threshold))
                    pred_pair_count = build_pair_counts(pred)

                    score = omega_index(
                        ground_truth_pair_counts=gt_pair_counts,
                        detected_pair_counts=pred_pair_count,
                        total_pairs=total_pairs,
                    )
                    sum_scores[s] += score

                    rewire_step(G=G, edge_stack=edge_stack, k=k_step, rng=rng)

                    dt_step = time.perf_counter() - t_step0
                    print(f"      [STEP {s + 1}/{steps}] dt={dt_step:.3f}s")

                # final point p=1.0
                pred = get_communities(G, closed_neighborhood_jaccard_predicate(threshold))
                pred_pair_count = build_pair_counts(pred)
                score = omega_index(
                    ground_truth_pair_counts=gt_pair_counts,
                    detected_pair_counts=pred_pair_count,
                    total_pairs=total_pairs,
                )
                sum_scores[steps] += score

                dt_run = time.perf_counter() - t_run0
                print(f"    [RUN {run_i}/{avg_times}] done  dt={dt_run:.2f}s")

            ys = [v / avg_times for v in sum_scores]

            # ---- Plot per-threshold ----
            plt.figure(figsize=(12, 6))
            plt.plot(xs, ys)

            plt.xlabel("p (fraction rewired)")
            plt.ylabel("Recognition success")

            plt.ylim(-0.05, 1.05)
            ticks01 = [i / 10 for i in range(11)]
            plt.xticks(ticks01)
            plt.yticks(ticks01)

            plt.grid(True, alpha=0.25)

            title = f"Rewiring vs Truth | n={n}, z={z}, rings={rings}, t={threshold:.3f}"
            plt.title(title)

            fname = (
                f"rewiring_vs_truth_n{n}_z{z}_r{rings}"
                f"_t{safe_sheet_num(threshold)}"
                f"_pstep{safe_sheet_num(p_step)}_avg{avg_times}.png"
            )
            path = os.path.join(out_dir, fname)

            plt.tight_layout()
            plt.savefig(path, dpi=220)
            plt.close()

            dt_thr = time.perf_counter() - t_thr0
            print(f"  [THR {ti}/{total_thresholds}] saved -> {path}  dt={dt_thr:.2f}s")

        dt_cfg = time.perf_counter() - t_cfg0
        print(f"[CONFIG {ci}/{total_configs}] done  dt={dt_cfg:.2f}s")

    dt_all = time.perf_counter() - t_global0
    print(f"\n[ALL DONE] total dt={dt_all:.2f}s")