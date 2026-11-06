"""
Controlled re-timing of the parameter sweeps on the full SNAP Amazon and DBLP graphs
one idle machine, one thread per run (OMP/MKL/OpenBLAS threads = 1),
every method in its own process reading the same edge list, the settings of big_job.py:

  hc        hc_scale.ScaleScores (enumeration and scoring of all wedges), then the cover at
            each of the 30 thresholds; the overlapping modularity of the 30 covers timed
            separately
  leiden    leidenalg, RB modularity, the 19 resolutions of big_job.LEIDEN_RB, one run each,
            n_iterations = -1 (to convergence), seed 1
  infomap   Infomap 2.10, two-level, one trial, the 13 Markov times of big_job.MARKOV
  lpa       networkx label propagation (parameter-free), one run

Reported per (network, method): wall time and CPU time of the sweep (graph construction
excluded: the time starts once the method's own input structure exists), the time of
the overlapping modularity of every cover (hc only), and the peak resident memory of
the process. Writes results/retime_big.csv.

    python retime_big.py [amazon dblp]
"""
import csv
import gzip
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = os.path.join(HERE, "..", "analysis", "snap")
OUT = os.path.join(HERE, "results", "retime_big.csv")
ENV1 = {k: "1" for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")}


def load(net):
    import numpy as np
    e = np.loadtxt(gzip.open(os.path.join(SNAP, f"com-{net}.ungraph.txt.gz"), "rt"), comments="#", dtype=np.int64)
    e = e[e[:, 0] != e[:, 1]]
    e = np.unique(np.sort(e, axis=1), axis=0)
    ids = np.unique(e)
    return np.searchsorted(ids, e), len(ids)


def peak_gb():
    import psutil
    mi = psutil.Process().memory_info()
    return getattr(mi, "peak_wset", getattr(mi, "rss", 0)) / 1e9


def child(net, method):
    import numpy as np
    sys.path.insert(0, HERE)
    from big_job import LEIDEN_RB, MARKOV, T_GRID
    e, n = load(net)
    res = dict(network=net, method=method)
    if method == "hc":
        import scipy.sparse as sp
        from hc_scale import ScaleScores
        from self_tune import overlapping_modularity
        A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n))
        w0, c0 = time.perf_counter(), time.process_time()
        S = ScaleScores(edges=e, n=n)
        S.nodes = list(range(n))
        build_w = time.perf_counter() - w0
        covers = [S.communities(t, variant="hc") for t in T_GRID]
        res.update(wall_s=time.perf_counter() - w0, cpu_s=time.process_time() - c0, build_s=build_w)
        w1 = time.perf_counter()
        idx = {i: i for i in range(n)}
        for cov in covers:
            overlapping_modularity(A, idx, cov)
        res["eq_s"] = time.perf_counter() - w1
    elif method == "leiden":
        import igraph as ig
        import leidenalg as la
        g = ig.Graph(n=n, edges=e.tolist())
        w0, c0 = time.perf_counter(), time.process_time()
        for gamma in LEIDEN_RB:
            la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter=gamma,
                              n_iterations=-1, seed=1)
        res.update(wall_s=time.perf_counter() - w0, cpu_s=time.process_time() - c0, runs=len(LEIDEN_RB))
    elif method == "infomap":
        from infomap import Infomap
        w0, c0 = time.perf_counter(), time.process_time()
        for mt in MARKOV:
            im = Infomap(f"--two-level --silent --seed 1 --num-trials 1 --markov-time {mt}")
            for a, b in e:
                im.add_link(int(a), int(b))
            im.run()
        res.update(wall_s=time.perf_counter() - w0, cpu_s=time.process_time() - c0, runs=len(MARKOV))
    elif method == "lpa":
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(map(tuple, e.tolist()))
        w0, c0 = time.perf_counter(), time.process_time()
        list(nx.community.label_propagation_communities(G))
        res.update(wall_s=time.perf_counter() - w0, cpu_s=time.process_time() - c0, runs=1)
    res["peak_gb"] = peak_gb()
    print("RESULT " + json.dumps(res), flush=True)


def main(nets):
    rows = []
    for net in nets:
        for method in ("hc", "lpa", "infomap", "leiden"):
            env = dict(os.environ, **ENV1)
            p = subprocess.run([sys.executable, "-u", __file__, "--child", net, method], env=env,
                               capture_output=True, text=True)
            line = [x for x in p.stdout.splitlines() if x.startswith("RESULT ")]
            if not line:
                print(net, method, "FAILED", p.stderr[-500:], flush=True)
                continue
            r = json.loads(line[0][7:])
            print(r, flush=True)
            rows.append(r)
            keys = sorted({k for x in rows for k in x})
            with open(OUT, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=keys)
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        child(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1:] or ["amazon", "dblp"])
