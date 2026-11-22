"""
BigClam and ego-splitting (bigclam.py, egosplit.py) on SNAP Amazon, DBLP and YouTube, full
graphs and the subgraphs induced by the ground-truth nodes: covers in big_job.py's format,
scored afterwards with the other methods by score_big.py. Seed 1, as big_job.py's Leiden and
Infomap. One task per (network, setting, K) for BigClam (K of methods2.bigclam_grid) and per (network, setting, local clustering) for ego-splitting (global
clustering: label propagation and Leiden-RB over its 13 resolutions, one run each).

    python run_more_big.py [--workers 4] [--only amazon,dblp] [--settings full,induced]
                           [--methods bigclam,egosplit] [--out results/big_dl/results_big]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SEED = 1


def _csr(e, n):
    import numpy as np
    A = __import__("big_job").adjacency(e, n)
    return A, A.indptr.astype(np.int64), A.indices.astype(np.int64)


def task(args):
    net, setting, method, prm, out = args
    from big_job import Store, load, to_arrays
    st = Store(out, net, setting, method, f"_{prm}")
    if st.done():
        return f"{time.strftime('%H:%M:%S')} skip  {net} {setting} {method} {prm}"
    t0 = time.time()
    try:
        e, ids = load(net, setting)
        A, indptr, indices = _csr(e, len(ids))
        if method == "bigclam":
            from bigclam import BigClam
            model = BigClam(indptr, indices, None, SEED)
            K = int(prm)
            import numpy as np
            model.rng = np.random.default_rng(SEED * 1009 + K)
            model.init(K)
            model.fit(seed=(SEED * 1009 + K) % (2**31))
            st.add(K, *to_arrays(model.communities()), A, ids, time.time() - t0)
        else:
            import networkx as nx
            from egosplit import egosplit_covers
            from methods2 import GRIDS2
            G = nx.Graph()
            G.add_nodes_from(range(len(ids)))
            G.add_edges_from(e.tolist())
            t1 = time.time()
            for param, cover in egosplit_covers(G, [g for g in GRIDS2["egosplit"] if g[0] == prm], SEED, restarts=1):
                st.add(param, *to_arrays(cover), A, ids, time.time() - t1)
                t1 = time.time()
        st.save()
        return f"{time.strftime('%H:%M:%S')} done  {net} {setting} {method} {prm} ({(time.time() - t0) / 60:.1f} min)"
    except Exception:
        msg = f"{net} {setting} {method} {prm}\n{traceback.format_exc()}"
    with open(os.path.join(out, "errors.txt"), "a") as fh:
        fh.write(f"more {msg}\n")
    return f"{time.strftime('%H:%M:%S')} ERROR {msg.splitlines()[0]} ({(time.time() - t0) / 60:.1f} min)"


def main():
    from big_job import load
    from methods2 import GRIDS2, bigclam_grid
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", default="amazon,dblp")
    ap.add_argument("--settings", default="full,induced")
    ap.add_argument("--methods", default="bigclam,egosplit")
    ap.add_argument("--out", default=os.path.join(HERE, "results", "big_dl", "results_big"))
    a = ap.parse_args()
    tasks = []
    for net in a.only.split(","):
        for s in a.settings.split(","):
            n = len(load(net, s)[1])
            for m in a.methods.split(","):
                prms = bigclam_grid(n) if m == "bigclam" \
                    else sorted({g[0] for g in GRIDS2["egosplit"]})
                tasks += [(net, s, m, p, a.out) for p in prms]
    tasks.sort(key=lambda x: (x[1] != "full", -(int(x[3]) if x[2] == "bigclam" else 1000)))   # the longest first
    print(f"{len(tasks)} tasks, {a.workers} workers", flush=True)
    with Pool(a.workers, maxtasksperchild=1) as pool:
        for i, m in enumerate(pool.imap_unordered(task, tasks), 1):
            print(f"[{i}/{len(tasks)}] {m}", flush=True)
    print(f"{time.strftime('%H:%M:%S')} ALL DONE", flush=True)


if __name__ == "__main__":
    main()
