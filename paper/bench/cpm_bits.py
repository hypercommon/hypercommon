"""
Clique percolation (Palla et al. 2005) for dense graphs with millions of maximal cliques,
with the output of networkx's k_clique_communities. Every maximal clique is a bit string
over the nodes; every pair of cliques is compared (AND and population count, compiled
with numba), and pairs sharing at least k-1 nodes are merged in a union-find, for every k
at once. Each thread keeps its own union-find over its rows; the partitions are merged at
the end, which gives the components of the union of all percolating pairs. Memory is
the bit strings plus one union-find per thread and k; time grows with the square of the
number of cliques (facebook_107: 2.2 million cliques).

Needs numba (a separate env; see paper/README.md):
    python cpm_bits.py <network> [k ...]       (real_more network; default k = 3..7)
    python cpm_bits.py --test
"""
from __future__ import annotations

import os
import sys
import time

import networkx as nx
import numba as nb
import numpy as np
from numba import types
from numba.extending import intrinsic


@intrinsic
def _popcount(typingctx, x):
    sig = types.uint64(types.uint64)

    def codegen(context, builder, signature, args):
        fn = builder.module.declare_intrinsic("llvm.ctpop", [args[0].type])
        return builder.call(fn, args)
    return sig, codegen


@nb.njit(cache=True, inline="always")
def _find(p, i):
    while p[i] != i:
        p[i] = p[p[i]]
        i = p[i]
    return i


@nb.njit(parallel=True, cache=True)
def _percolate(bits, sizes, ks, nthreads):
    m, W = bits.shape
    K = len(ks)
    parents = np.empty((nthreads, K, m), np.int32)
    for t in range(nthreads):
        for kk in range(K):
            for i in range(m):
                parents[t, kk, i] = i
    kmin = ks[0]
    for t in nb.prange(nthreads):
        P = parents[t]
        for i in range(t, m, nthreads):            # interleaved rows balance the triangle
            si = sizes[i]
            for j in range(i + 1, m):
                ov = 0
                for w in range(W):
                    ov += _popcount(bits[i, w] & bits[j, w])
                if ov < kmin - 1:
                    continue
                sj = sizes[j]
                for kk in range(K):
                    k = ks[kk]
                    if ov >= k - 1 and si >= k and sj >= k:
                        a = _find(P[kk], i)
                        b = _find(P[kk], j)
                        if a != b:
                            if a < b:
                                P[kk, b] = a
                            else:
                                P[kk, a] = b
    # merge the threads' partitions
    out = np.empty((K, m), np.int32)
    for kk in range(K):
        g = out[kk]
        for i in range(m):
            g[i] = i
        for t in range(nthreads):
            for i in range(m):
                r = _find(parents[t, kk], i)
                if r != i:
                    a = _find(g, i)
                    b = _find(g, r)
                    if a != b:
                        if a < b:
                            g[b] = a
                        else:
                            g[a] = b
        for i in range(m):
            g[i] = _find(g, i)
    return out


def k_clique_covers(G, ks=(3, 4, 5, 6, 7), nthreads=None, log=False):
    ks = sorted(ks)
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    t0 = time.time()
    cl = [np.fromiter((idx[u] for u in c), np.int64, len(c)) for c in nx.find_cliques(G) if len(c) >= ks[0]]
    m = len(cl)
    out = {k: [] for k in ks}
    if m == 0:
        return out
    sizes = np.array([len(c) for c in cl], np.int32)
    W = (n + 63) // 64
    bits = np.zeros((m, W), np.uint64)
    rows = np.repeat(np.arange(m), sizes)
    cols = np.concatenate(cl)
    np.bitwise_or.at(bits, (rows, cols // 64), np.left_shift(np.uint64(1), (cols % 64).astype(np.uint64)))
    if log:
        print(f"{m} maximal cliques of >= {ks[0]} nodes, {W} words each ({time.time() - t0:.0f} s)", flush=True)
    nthreads = nthreads or nb.get_num_threads()
    t1 = time.time()
    lab = _percolate(bits, sizes, np.array(ks, np.int64), nthreads)
    if log:
        print(f"percolation of {m * (m - 1) / 2:.3g} pairs on {nthreads} threads: {time.time() - t1:.0f} s", flush=True)
    for kk, k in enumerate(ks):
        ok = np.flatnonzero(sizes >= k)
        if not len(ok):
            continue
        rr = np.repeat(lab[kk][ok], sizes[ok])
        cc = np.concatenate([cl[i] for i in ok])
        pairs = np.unique(np.c_[rr, cc], axis=0)
        split = np.flatnonzero(np.diff(pairs[:, 0])) + 1
        out[k] = [frozenset(nodes[j] for j in grp) for grp in np.split(pairs[:, 1], split)]
    return out


def _test():
    import random
    for s in range(40):
        rnd = random.Random(s)
        G = nx.gnp_random_graph(rnd.randint(20, 200), rnd.uniform(0.05, 0.35), seed=s)
        if s % 2:
            G = nx.relabel_nodes(G, {u: f"v{u}" for u in G})
        got = k_clique_covers(G, (3, 4, 5, 6, 7), nthreads=1 + s % 5)
        for k, cov in got.items():
            ref = {frozenset(c) for c in nx.community.k_clique_communities(G, k)}
            assert set(cov) == ref and len(cov) == len(ref), (s, k)
    print("cpm_bits matches networkx on 40 random graphs")


if __name__ == "__main__":
    if sys.argv[1:] == ["--test"]:
        _test()
        sys.exit()
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(1, os.path.join(here, "..", "analysis"))
    from real_more import load
    args = sys.argv[1:]
    save = None
    if "--save" in args:
        i = args.index("--save")
        save = args[i + 1]
        del args[i:i + 2]
    G, _ = load(args[0])
    ks = [int(x) for x in args[1:]] or [3, 4, 5, 6, 7]
    t0 = time.time()
    cov = k_clique_covers(G, ks, log=True)
    secs = time.time() - t0
    print({k: len(v) for k, v in cov.items()}, f"{secs:.0f} s", flush=True)
    if save:
        import pickle
        with open(save, "wb") as f:
            pickle.dump({"covers": {k: [sorted(c) for c in v] for k, v in cov.items()}, "secs": secs}, f)
