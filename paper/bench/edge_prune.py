"""Wedges left to score at threshold t: all, after degree-ratio pruning, after exact edge-J pruning
(the numbers of the Discussion).

    python bench/edge_prune.py <edge list> ...      (SNAP .txt.gz, or a graph folder's edges.tsv)
Prints a table per graph and writes figures/data/edge_prune.csv (one row per graph and t).
"""
import gzip
import os
import sys

import numpy as np
import scipy.sparse as sp


def load(path):
    if path.endswith(".gz"):
        rows = [l.split() for l in gzip.open(path, "rt") if not l.startswith("#")]
    else:
        rows = [l.split()[:2] for l in open(path) if l.strip() and not l.startswith(("#", "source", "u\t"))]
    u = np.array([r[0] for r in rows]); v = np.array([r[1] for r in rows])
    ids, inv = np.unique(np.r_[u, v], return_inverse=True)
    a, b = inv[:len(u)], inv[len(u):]
    n = len(ids)
    A = sp.csr_matrix((np.ones(2 * len(a), np.int8), (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    A.sum_duplicates(); A.data[:] = 1; A.setdiag(0); A.eliminate_zeros()
    return A.tocsr()


ROWS = []


def main(path):
    A = load(path)
    d = np.asarray(A.sum(1)).ravel().astype(np.int64)
    U = sp.triu(A, 1).tocoo()
    u, v = U.row, U.col
    # common neighbours per edge = (A @ A)[u, v], row blocks to save memory
    c = np.zeros(len(u), np.int64)
    order = np.argsort(u, kind="stable")
    Ai = A.astype(np.int32)
    B = 20000
    for s in range(0, A.shape[0], B):
        P = (Ai[s:s + B] @ Ai).tocsr()
        sel = order[(u[order] >= s) & (u[order] < s + B)]
        c[sel] = np.asarray(P[u[sel] - s, v[sel]]).ravel()
    J = (c + 2) / (d[u] + d[v] - c)
    ratio = (np.minimum(d[u], d[v]) + 1) / (np.maximum(d[u], d[v]) + 1)
    W = int((d * (d - 1) // 2).sum())
    name = path.replace(os.sep, "/").split("/")[-2] if "edges" in path else os.path.basename(path)
    print(f"{name}: n={A.shape[0]} m={len(u)} wedges={W:,}")
    print("    t   ratio-pruned   edge-J-pruned   edges kept (J)")
    for t in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
        out = []
        for keep in (ratio >= t, J >= t):
            dk = np.bincount(np.r_[u[keep], v[keep]], minlength=A.shape[0]).astype(np.int64)
            out.append(int((dk * (dk - 1) // 2).sum()) / W)
        print(f"  {t:.1f}   {out[0]:12.1%}   {out[1]:13.1%}   {(J >= t).mean():8.1%}")
        ROWS.append(dict(graph=name, n=A.shape[0], m=len(u), wedges=W, t=t, ratio_pruned=out[0],
                         edge_pruned=out[1], edges_kept=float((J >= t).mean())))


if __name__ == "__main__":
    import pandas as pd
    for p in sys.argv[1:]:
        main(p)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "data", "edge_prune.csv")
    pd.DataFrame(ROWS).to_csv(out, index=False)
