"""
HyperCommon on the full SNAP YouTube graph, for a machine with ~128 GB of memory.

Builds hc_big.BigScores with degree-ratio pruning at t_min (exact for every t >= t_min)
and, for every variant and every t of the grid, writes the cover and its overlapping
modularity. Covers are scored against the ground truth later, at home.

    python youtube_job.py [--tmin 0.02] [--max-edges N] [--out DIR]

Output in DIR (default: results next to this file):
  summary.csv            one row per (variant, t): eq, communities, memberships, seconds
  covers_<variant>.npz   per t: t<t>_node (SNAP ids, grouped by community), t<t>_size
  meta.json              graph and build statistics, peak memory, versions
  log.txt                progress
Reruns skip variants whose covers file exists.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import gzip
import json
import os
import platform
import sys
import time
import warnings

import numpy as np
import scipy
import scipy.sparse as sp

warnings.filterwarnings("ignore", message="Input line")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hc_big import BigScores, overlapping_modularity  # noqa: E402

T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]
VARIANTS = ["hc", "hc_tri", "hc_2test", "hc_econn"]


def peak_gb():
    if os.name == "nt":
        class PMC(ctypes.Structure):
            _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + \
                       [(f, ctypes.c_size_t) for f in ("PeakWorkingSetSize", "WorkingSetSize",
                                                       "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage",
                                                       "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage",
                                                       "PagefileUsage", "PeakPagefileUsage")]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.windll.kernel32
        k32.GetCurrentProcess.restype = ctypes.c_void_p
        k32.K32GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(PMC), ctypes.c_ulong]
        k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return round(pmc.PeakWorkingSetSize / 2**30, 2)
    import resource
    # ru_maxrss is in bytes on macOS, in kilobytes on Linux
    unit = 2**30 if sys.platform == "darwin" else 2**20
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / unit, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default=os.path.join(HERE, "data", "com-youtube.ungraph.txt.gz"))
    ap.add_argument("--tmin", type=float, default=0.02)
    ap.add_argument("--max-edges", type=int, default=0, help="test run on the first N edges")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")
    t_start = time.time()

    def log(msg):
        line = f"{time.strftime('%H:%M:%S')} [{(time.time() - t_start) / 60:6.1f} min, peak {peak_gb()} GB] {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    log(f"start: python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
        f"{os.cpu_count()} cpus, t_min {a.tmin}, max_edges {a.max_edges or 'all'}")
    with gzip.open(a.edges, "rt") as fh:
        raw = np.loadtxt(fh, comments="#", dtype=np.int64, max_rows=a.max_edges or None)
    ids, e = np.unique(raw, return_inverse=True)        # SNAP ids -> 0..n-1
    e = e.reshape(-1, 2)
    n = len(ids)
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n))
    A.sum_duplicates(); A.data[:] = 1; A.setdiag(0); A.eliminate_zeros()
    log(f"graph: {n} nodes, {A.nnz // 2} edges")

    t0 = time.time()
    B = BigScores(e, n, t_min=a.tmin, log=log)
    build = round(time.time() - t0, 1)
    meta = dict(nodes=n, edges=A.nnz // 2, t_min=a.tmin, max_edges=a.max_edges, wedges_full=B.n_wedges_full,
                wedges_kept=B.n_wedges, triples=len(B.triples), pairs=B.n_pairs, build_seconds=build,
                peak_gb_after_build=peak_gb(), python=platform.python_version(), numpy=np.__version__,
                scipy=scipy.__version__, cpus=os.cpu_count(), machine=platform.platform())
    log(f"built in {build} s")

    sfile = os.path.join(a.out, "summary.csv")
    new = not os.path.exists(sfile)
    sf = open(sfile, "a", newline="")
    w = csv.DictWriter(sf, fieldnames=["variant", "t", "eq", "n_comms", "memberships", "nodes_covered", "seconds"])
    if new:
        w.writeheader()
    for v in VARIANTS:
        cfile = os.path.join(a.out, f"covers_{v}.npz")
        if os.path.exists(cfile):
            log(f"{v}: done before, skipped")
            continue
        arrs = {}
        for t in T_GRID:
            if t < a.tmin:
                continue
            t1 = time.time()
            node, comm = B.cover(t, v)
            eq = overlapping_modularity(A, node, comm)
            secs = round(time.time() - t1, 1)
            size = np.bincount(comm) if len(comm) else np.zeros(0, np.int64)
            arrs[f"t{t:.2f}_node"] = ids[node].astype(np.int32)
            arrs[f"t{t:.2f}_size"] = size.astype(np.int32)
            w.writerow(dict(variant=v, t=t, eq=round(eq, 6), n_comms=len(size), memberships=len(node),
                            nodes_covered=len(np.unique(node)), seconds=secs))
            sf.flush()
            log(f"{v} t={t:.2f}: {len(size)} communities, eq {eq:.4f}, {secs} s")
        np.savez_compressed(cfile + ".tmp.npz", **arrs)
        os.replace(cfile + ".tmp.npz", cfile)
    meta["peak_gb"] = peak_gb()
    meta["total_minutes"] = round((time.time() - t_start) / 60, 1)
    with open(os.path.join(a.out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    log("ALL DONE. Please send back the whole results folder.")


if __name__ == "__main__":
    main()
