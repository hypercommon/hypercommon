"""
Running time against network size, single process, on an otherwise idle machine.

Graphs: overlapping ring lattices (rings of 100, k = 8 and 16, o_n = 0.05, p = 0.2)
and overlapping LFR graphs (<k> = 20, k_max = 50, sizes 20-100, o_n = 0.1, mu = 0.3),
n from 1e3 to 2e5 (LFR to 1e5). Three repetitions each.

Timed: HyperCommon build (wedge enumeration + all triple scores, hc_scale), one
threshold's communities, the label-free threshold sweep (30 thresholds + overlapping
modularity), a LOCAL query (the community of one random seed node, BFS over shared
pairs with on-demand commonality, hypercommon.algorithm.get_node_community), Leiden
(gamma = 1) and Infomap. W (wedges) and peak traced memory of the HyperCommon build
are recorded. Output: results/runtime.csv
"""
import csv
import os
import random
import sys
import time
import tracemalloc

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graphs                                                     # noqa: E402
from hc_scale import ScaleScores                                  # noqa: E402
from self_tune import overlapping_modularity                      # noqa: E402
from hypercommon.algorithm import get_node_community             # noqa: E402  (repo root via graphs)
from predicates.jaccard import closed_neighborhood_jaccard_predicate  # noqa: E402

OUT = os.path.join(HERE, "results", "runtime.csv")
T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]


def ring_graph(n, k, seed):
    r = n // 100
    return graphs.rings([100] * r, [k] * r, 0.05, 0.2, seed)[0]


def lfr_graph(n, seed):
    return graphs.lfr(dict(N=n, k=20, maxk=50, t1=2, t2=1, minc=20, maxc=100, on_frac=0.1, om=2), 0.3, seed)[0]


def timed(f):
    t = time.perf_counter()
    out = f()
    return out, time.perf_counter() - t


def main():
    cases = [("rings_k8", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)]
    cases += [("rings_k16", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)]
    cases += [("lfr", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000)]
    new = not os.path.exists(OUT)
    done = set()
    if not new:     # resumable: skip (family, n, rep) already recorded
        with open(OUT, newline="") as rf:
            done = {(r["family"], int(r["n"]), int(r["rep"])) for r in csv.DictReader(rf)}
    fh = open(OUT, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=["family", "n", "m", "W", "rep", "hc_build_s", "hc_one_t_s",
                                       "hc_sweep_s", "hc_build_peak_gb", "hc_local_s", "leiden_s", "infomap_s"])
    if new:
        w.writeheader()
    for fam, n in cases:
        for rep in range(3):
            G = ring_graph(n, 8 if fam == "rings_k8" else 16, rep) if fam.startswith("rings") else lfr_graph(n, rep)
            if (fam, G.number_of_nodes(), rep) in done:
                continue
            nodes = list(G.nodes())
            idx = {u: i for i, u in enumerate(nodes)}
            e = np.array([(idx[a], idx[b]) for a, b in G.edges()], dtype=np.int64)
            # timing without the memory tracer (as for the competitors) ...
            S, t_build = timed(lambda: ScaleScores(edges=e, n=len(nodes)))
            # ... and peak memory in a separate, untimed build, on the first repetition only
            peak = float("nan")
            if rep == 0:
                del S
                tracemalloc.start()
                S = ScaleScores(edges=e, n=len(nodes))
                peak = tracemalloc.get_traced_memory()[1] / 1e9
                tracemalloc.stop()
            _, t_one = timed(lambda: S.communities(0.12))
            A = S.A.astype(float)
            ident = {i: i for i in range(len(nodes))}

            def sweep():
                best = (-np.inf, None)
                for t in T_GRID:
                    q = overlapping_modularity(A, ident, S.communities(t))
                    best = max(best, (q, t))
                return best
            _, t_sweep = timed(sweep)
            seed_node = random.Random(rep).choice(nodes)
            _, t_local = timed(lambda: get_node_community(G, closed_neighborhood_jaccard_predicate(0.12), seed_node))
            import igraph as ig
            import leidenalg as la
            g = ig.Graph(n=len(nodes), edges=e.tolist())
            _, t_leiden = timed(lambda: la.find_partition(g, la.RBConfigurationVertexPartition,
                                                          resolution_parameter=1.0, seed=1))
            from infomap import Infomap

            def run_infomap():
                im = Infomap("--two-level --silent --seed 1")
                for a, b in e.tolist():
                    im.add_link(a, b)
                im.run()
            _, t_inf = timed(run_infomap)
            row = dict(family=fam, n=len(nodes), m=len(e), W=S.n_wedges, rep=rep, hc_build_s=round(t_build, 3),
                       hc_one_t_s=round(t_one, 3), hc_sweep_s=round(t_sweep, 3), hc_build_peak_gb=round(peak, 3),
                       hc_local_s=round(t_local, 4), leiden_s=round(t_leiden, 3), infomap_s=round(t_inf, 3))
            w.writerow(row)
            fh.flush()
            print(row, flush=True)


if __name__ == "__main__":
    main()
