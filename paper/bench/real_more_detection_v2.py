"""
Detection on the additional real networks (analysis/real_more.py): ten Facebook
ego networks with overlapping social circles, political blogs, a primary school,
college football and political books.

Every method of the benchmark runs its whole grid (original + extensions; SLPA
with 100 iterations; Leiden over a wide resolution range) and HyperCommon runs with
every commonality function of the family; every HyperCommon variant searches the same
100 quantile levels of its own triple scores. Stochastic methods run with seeds 0..4.
Evaluation is on the nodes that carry at least one recorded group, to which every
detected community is restricted (as for the SNAP networks); every nonempty restricted
community is kept, and the *_ge3 columns also drop those of fewer than three nodes.

Per (network, method, param, seed): F1 (best-match), ONMI, Omega, coverage.
Label-free rows: HyperCommon (each f) at the threshold maximising overlapping
modularity; other methods at their library defaults.

    python real_more_detection_v2.py [workers]
Output: results/real_more_v2.csv (appended, resumable by (network, method)).
"""
from __future__ import annotations

import csv
import os
import sys
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "analysis"))

OUT = os.path.join(HERE, "results", "real_more_v2.csv")
# Two scoring rules, both recorded: "all" keeps every nonempty detected community after
# restriction to the labelled nodes (primary); "ge3" also drops restricted communities
# with fewer than three nodes (sensitivity check; ground-truth groups have >= 3).
FIELDS = ["network", "method", "param", "seed", "label_free",
          "f1", "onmi", "omega", "coverage", "n_comms",
          "f1_ge3", "onmi_ge3", "omega_ge3", "coverage_ge3", "n_comms_ge3",
          "seconds", "n_nodes", "n_edges", "error"]
# Every HyperCommon variant and every commonality function searches the same grid:
# 100 quantile levels of its own triple scores (q = 0 admits every wedge).
Q100 = [round(0.01 * i, 2) for i in range(100)]
LEIDEN_GRID = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]
DEFAULTS = {"leiden": 1.0, "slpa": 0.1, "demon": 0.25, "angel": 0.25, "cpm": 3, "scan": (0.5, 2), "scan_j": (0.3, 2)}
STOCHASTIC = {"leiden", "infomap", "lpa", "slpa", "demon", "angel"}
METHODS = (["components", "hc", "hc_tri", "hc_2test"]
           + ["hcf_" + f for f in ("jaccard", "cosine", "max_norm", "min_norm", "i2_over_u", "count", "resource_alloc")]
           + ["cpm", "leiden", "infomap", "lpa", "walktrap", "slpa", "demon", "angel", "scan", "scan_j", "linkcomm"])


# clique percolation is out of reach on facebook_1912 (as real_more_detection_v3.SKIP); its
# facebook_107 rows come from cpm_bits.py through cpm_rows_from_bits.py
SKIP = {("facebook_107", "cpm"), ("facebook_1912", "cpm")}


def task(args):
    net, method = args
    import networkx as nx
    import numpy as np
    import cover_metrics as M
    import methods as ME
    from real_more import load
    from self_tune import overlapping_modularity, _adj

    G, truth = load(net)
    universe = set().union(*truth)
    nodes_u = list(universe)

    Gu = G.subgraph(universe)

    def restrict(cov, min_size):
        c = [set(x) & universe for x in cov]
        return [x for x in c if len(x) >= min_size]

    rows = []

    def scores(r, suffix):
        return {f"f1{suffix}": M.f1_best_match(r, truth, universe), f"onmi{suffix}": M.onmi(r, truth, Gu),
                f"omega{suffix}": M.omega(r, truth, nodes_u), f"coverage{suffix}": M.coverage(r, truth),
                f"n_comms{suffix}": len(r)}

    def emit(param, cov, secs, seed=0, label_free=""):
        rows.append(dict(network=net, method=method, param=param, seed=seed, label_free=label_free,
                         **scores(restrict(cov, 1), ""), **scores(restrict(cov, 3), "_ge3"),
                         seconds=round(secs, 3), n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(), error=""))

    try:
        t0 = time.time()
        if method == "components":
            emit("", [set(c) for c in nx.connected_components(G)], time.time() - t0, label_free="trivial")
        elif method.startswith("hc"):
            grid = Q100
            if method.startswith("hcf_"):
                from hc_family import FamilyScores
                ts = FamilyScores(G)
                f = method[4:]
                get = lambda q: ts.communities_f(f, q)
            else:
                from hc_fast import TripleScores
                ts = TripleScores(G)
                score = {"hc": ts.s3, "hc_2test": ts.s2, "hc_tri": ts.s3[ts.closed]}[method]

                def get(q, ts=ts, score=score):
                    thr = np.quantile(score, q) if (q > 0 and len(score)) else -np.inf
                    return ts.communities(method, thr)
            A, idx = _adj(G)
            best = (-np.inf, None, None)
            for t in grid:
                t1 = time.time()
                cov = get(t)
                emit(t, cov, time.time() - t1)
                q = overlapping_modularity(A, idx, cov)
                if q > best[0]:
                    best = (q, t, cov)
            emit(best[1], best[2], 0.0, label_free="max_overlapping_modularity")
        else:
            seeds = range(5) if method in STOCHASTIC else [0]
            for seed in seeds:
                if method == "leiden":
                    import igraph as ig
                    import leidenalg as la
                    nl = list(G.nodes())
                    ix = {u: i for i, u in enumerate(nl)}
                    g = ig.Graph(n=len(nl), edges=[(ix[a], ix[b]) for a, b in G.edges()])
                    for gamma in LEIDEN_GRID:
                        t1 = time.time()
                        part = la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter=gamma,
                                                 seed=seed)
                        emit(gamma, [{nl[i] for i in c} for c in part], time.time() - t1, seed,
                             "default" if gamma == DEFAULTS["leiden"] else "")
                    continue
                name = "slpa100" if method == "slpa" else method
                for param, cov in ME.run_method(name, G, None, seed):
                    emit(param, cov, 0.0, seed, "default" if (param == DEFAULTS.get(method) or param is None) else "")
                if method in ME.GRIDS_X:
                    for param, cov in ME.run_method(method + "_x", G, None, seed):
                        emit(param, cov, 0.0, seed)
    except Exception as e:
        rows.append(dict(network=net, method=method, param="", seed="", label_free="", f1="", onmi="", omega="",
                         coverage="", n_comms="", seconds="", n_nodes="", n_edges="",
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=3)!r}"[:600]))
    return rows


def main():  # noqa: C901
    from real_more import ALL
    nets = [n for n in ALL if n != "cora"]
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            done = {(r["network"], r["method"]) for r in csv.DictReader(fh) if not r["error"]}
    tasks = [(n, m) for n in nets for m in METHODS if (n, m) not in done and (n, m) not in SKIP]
    # big, slow tasks first
    order = {"linkcomm": 0, "slpa": 1, "demon": 2, "angel": 3, "walktrap": 4}
    tasks.sort(key=lambda t: order.get(t[1], 9))
    print(f"{len(tasks)} tasks", flush=True)
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as fh, Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        for i, rows in enumerate(pool.imap_unordered(task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            print(f"  {i}/{len(tasks)} {rows[0]['network']} {rows[0]['method']} {rows[0]['error'][:120]}", flush=True)


if __name__ == "__main__":
    main()
