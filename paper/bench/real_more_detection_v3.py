"""
Detection on the fourteen Netzschleuder networks, main runs.

Changes from real_more_detection_v2.py: the method set and grids of methods2.py
(Leiden with RB modularity and with the Constant Potts Model, Infomap over Markov
times, link communities over their similarity cut, the ablations hc_2test and
hc_econn), and every cover also scored with `eq`, its overlapping modularity on the
whole graph, so that every method's label-free choice -- its parameter with the
largest eq -- is read from the same rows. HyperCommon variants and link communities
search 100 quantile levels of their own scores; HC with each commonality function as
before. OSLOM (the OSLOM2 binary, run_oslom.py; its threshold t over T_OSLOM, five seeds) was
added after the first run, which a rerun picks up as the only missing method. Scoring
as in v2: detected communities restricted to the labelled nodes, every
nonempty one kept (f1, ...) or only those of at least three nodes (*_ge3).

    python real_more_detection_v3.py [workers]
Output: results/real_more_v3.csv (appended, resumable by (network, method)).
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

OUT = os.path.join(HERE, "results", "real_more_v3.csv")
FIELDS = ["network", "method", "param", "seed", "eq",
          "f1", "onmi", "omega", "coverage", "n_comms",
          "f1_ge3", "onmi_ge3", "omega_ge3", "coverage_ge3", "n_comms_ge3",
          "seconds", "n_nodes", "n_edges", "error"]
Q100 = [round(0.01 * i, 2) for i in range(100)]
STOCHASTIC = {"leiden", "leiden_cpm", "infomap", "lpa", "slpa", "demon", "angel", "oslom"}
FUNCS = ("jaccard", "cosine", "max_norm", "min_norm", "i2_over_u", "count", "resource_alloc")
METHODS = (["components", "hc", "hc_tri", "hc_2test", "hc_econn", "linkcomm"]
           + ["hcf_" + f for f in FUNCS]
           + ["cpm", "leiden", "leiden_cpm", "infomap", "lpa", "walktrap", "slpa", "demon", "angel",
              "scan", "scan_j", "oslom"])
# clique percolation: facebook_1912 has more than 14 million maximal cliques (up to 68 nodes),
# out of reach; facebook_107 has 2.2 million, too slow for bench/cpm_fast.py, and its rows come
# from cpm_bits.py (10.6 h on 12 threads) through cpm_rows_from_bits.py
SKIP = {("facebook_107", "cpm"), ("facebook_1912", "cpm")}


def task(args):
    net, method = args
    import networkx as nx
    import numpy as np
    import cover_metrics as M
    from methods2 import run_method2
    from real_more import load
    from self_tune import overlapping_modularity, _adj

    G, truth = load(net)
    universe = set().union(*truth)
    nodes_u = list(universe)
    Gu = G.subgraph(universe)
    A, idx = _adj(G)

    def restrict(cov, min_size):
        c = [set(x) & universe for x in cov]
        return [x for x in c if len(x) >= min_size]

    def scores(r, suffix):
        return {f"f1{suffix}": M.f1_best_match(r, truth, universe), f"onmi{suffix}": M.onmi(r, truth, Gu),
                f"omega{suffix}": M.omega(r, truth, nodes_u), f"coverage{suffix}": M.coverage(r, truth),
                f"n_comms{suffix}": len(r)}

    rows = []

    def emit(param, cov, secs, seed=0):
        rows.append(dict(network=net, method=method, param=param, seed=seed,
                         eq=overlapping_modularity(A, idx, cov),
                         **scores(restrict(cov, 1), ""), **scores(restrict(cov, 3), "_ge3"),
                         seconds=round(secs, 3), n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(),
                         error=""))

    try:
        t0 = time.time()
        if method == "components":
            emit("", [set(c) for c in nx.connected_components(G)], time.time() - t0)
        elif method.startswith("hcf_"):
            from hc_family import FamilyScores
            ts = FamilyScores(G)
            for q in Q100:
                t1 = time.time()
                emit(q, ts.communities_f(method[4:], q), time.time() - t1)
        elif method in ("hc", "hc_tri", "hc_2test", "hc_econn", "linkcomm"):
            from hc_fast import TripleScores
            ts = TripleScores(G)
            if method == "linkcomm":       # the link similarity of every linked edge pair
                score = np.concatenate([js[ok] for _, ok, js in ts._links])
            else:
                score = {"hc": ts.s3, "hc_econn": ts.s3, "hc_2test": ts.s2, "hc_tri": ts.s3[ts.closed]}[method]
            for q in Q100:
                thr = np.quantile(score, q) if (q > 0 and len(score)) else -np.inf
                t1 = time.time()
                cov = ts.linkcomm(thr) if method == "linkcomm" else ts.communities(method, thr)
                emit(q, cov, time.time() - t1)
        else:
            for seed in (range(5) if method in STOCHASTIC else [0]):
                t1 = time.time()
                for param, cov in run_method2(method, G, None, seed):
                    emit(param, cov, time.time() - t1, seed)
                    t1 = time.time()
    except Exception as e:
        rows.append(dict(network=net, method=method,
                         error=f"{type(e).__name__}: {e} | {traceback.format_exc(limit=3)!r}"[:600]))
    return rows


def main():
    from real_more import ALL
    nets = [n for n in ALL if n != "cora"]
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            done = {(r["network"], r["method"]) for r in csv.DictReader(fh) if not r["error"]}
    tasks = [(n, m) for n in nets for m in METHODS if (n, m) not in done and (n, m) not in SKIP]
    order = {"slpa": 0, "demon": 1, "angel": 2, "walktrap": 3, "cpm": 4}
    tasks.sort(key=lambda t: order.get(t[1], 9))
    print(f"{len(tasks)} tasks", flush=True)
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as fh, Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        for i, rows in enumerate(pool.imap_unordered(task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            print(f"  {i}/{len(tasks)} {rows[0]['network']} {rows[0]['method']} {rows[0].get('error', '')[:120]}",
                  flush=True)


if __name__ == "__main__":
    main()
