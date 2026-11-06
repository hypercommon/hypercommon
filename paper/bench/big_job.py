"""
Covers of every method on SNAP Amazon, DBLP and YouTube (full graphs and the subgraphs
induced by the ground-truth nodes), for a large machine. Scoring against the ground
truth is done afterwards (score_big.py); this job only produces covers.

Methods, as in real_detection_v4.py:
  HyperCommon hc, hc_tri, hc_2test, hc_econn and link communities on T_GRID (ScaleScores);
  HyperCommon with each commonality function of the family on 50 quantile levels;
  Leiden (RB and CPM), Infomap over Markov times, label propagation, connected components;
  on the induced subgraphs also clique percolation and SLPA.
On YouTube's full graph only the rivals run here: its HyperCommon covers come from
youtube_job.py (degree-ratio pruning, exact for Jaccard only).

    python big_job.py [--workers 16] [--only amazon,dblp,youtube] [--settings full,induced] [--out DIR]

Every task writes DIR/<net>_<setting>_<method>[_<param>].npz (per parameter p: p<p>_node
with SNAP ids grouped by community, p<p>_size) and a matching .csv with eq, community
count, seconds. Finished tasks are skipped on a rerun.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import time
import traceback
import warnings
from multiprocessing import Pool

import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")

T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]
Q50 = [round(0.02 * i, 2) for i in range(50)]
LEIDEN_RB = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
             256.0, 512.0, 1024.0, 2048.0, 4096.0]
LEIDEN_CPM = [float(f"{x:.3g}") for x in np.geomspace(1e-5, 0.5, 20)]
MARKOV = [float(f"{x:.4g}") for x in np.geomspace(0.25, 16, 13)]
FUNCS = ("jaccard", "cosine", "max_norm", "min_norm", "i2_over_u", "count", "resource_alloc")
FUNCTIONS = {                                   # as hc_family.FUNCTIONS
    "jaccard": lambda P: P["I"] / P["U"],
    "count": lambda P: P["I"].astype(float),
    "i2_over_u": lambda P: P["I"] ** 2 / P["U"],
    "cosine": lambda P: P["I"] / np.sqrt(P["su"] * P["sv"]),
    "min_norm": lambda P: P["I"] / np.minimum(P["su"], P["sv"]),
    "max_norm": lambda P: P["I"] / np.maximum(P["su"], P["sv"]),
    "resource_alloc": lambda P: P["ra"],
}


# ---------------------------------------------------------------------------- data
def load(net, setting):
    """Edges as int64 (k, 2) over compact ids, the SNAP id of each compact id."""
    e = np.loadtxt(gzip.open(os.path.join(DATA, f"com-{net}.ungraph.txt.gz"), "rt"), comments="#", dtype=np.int64)
    e = e[e[:, 0] != e[:, 1]]
    e = np.unique(np.sort(e, axis=1), axis=0)
    if setting == "induced":
        cfile = "all" if net == "youtube" else "top5000"
        uni = set()
        with gzip.open(os.path.join(DATA, f"com-{net}.{cfile}.cmty.txt.gz"), "rt") as fh:
            for line in fh:
                c = line.split()
                if len(c) >= 3:
                    uni.update(int(x) for x in c)
        uni = np.array(sorted(uni), np.int64)
        e = e[np.isin(e[:, 0], uni) & np.isin(e[:, 1], uni)]
        ids = uni                               # isolated universe nodes stay nodes
    else:
        ids = np.unique(e)
    comp = np.searchsorted(ids, e)
    return comp, ids


def adjacency(e, n):
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n))
    A.sum_duplicates()
    A.data[:] = 1
    return A


def to_arrays(cover):
    """List of sets of compact ids -> (node, comm) arrays."""
    cover = [c for c in cover if len(c)]
    if not cover:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    node = np.fromiter((u for c in cover for u in c), np.int64)
    comm = np.repeat(np.arange(len(cover)), [len(c) for c in cover])
    return node, comm


class Store:
    def __init__(self, out, net, setting, method, tag=""):
        self.base = os.path.join(out, f"{net}_{setting}_{method}{tag}")
        self.net, self.setting, self.method = net, setting, method
        self.arrs, self.rows = {}, []

    def done(self):
        return os.path.exists(self.base + ".npz")

    def add(self, param, node, comm, A, ids, secs, **extra):
        from hc_big import overlapping_modularity
        order = np.lexsort((node, comm))
        node, comm = node[order], comm[order]
        size = np.bincount(comm) if len(comm) else np.zeros(0, np.int64)
        self.arrs[f"p{param}_node"] = ids[node].astype(np.int64)
        self.arrs[f"p{param}_size"] = size.astype(np.int64)
        self.rows.append(dict(network=self.net, setting=self.setting, method=self.method, param=param,
                              eq=round(overlapping_modularity(A, node, comm), 6), n_comms=len(size),
                              memberships=len(node), seconds=round(secs, 2), **extra))

    def save(self):
        with open(self.base + ".csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(self.rows[0]) if self.rows else ["network"])
            w.writeheader()
            for r in self.rows:
                w.writerow(r)
        np.savez_compressed(self.base + ".tmp.npz", **self.arrs)
        os.replace(self.base + ".tmp.npz", self.base + ".npz")


# ---------------------------------------------------------------------------- tasks
def pair_quantities(S):
    """Raw neighbourhood quantities of every pair of S's triples (as hc_family): closed
    intersection I, union U, closed sizes su, sv, resource allocation ra."""
    A, n, d = S.A, S.n, S.deg
    pu, pv = S.pair_nodes[:, 0].astype(np.int64), S.pair_nodes[:, 1].astype(np.int64)
    keys = pu * n + pv
    order = np.argsort(keys)
    keys, pu_s = keys[order], pu[order]
    Aw = (A.multiply(1.0 / np.maximum(d, 1)[None, :])).tocsr().astype(np.float64)   # A D^-1
    c = np.zeros(len(keys))
    ra = np.zeros(len(keys))
    rows = np.unique(pu_s)
    cost = np.cumsum((A @ d.astype(np.float64))[rows])
    start = 0
    while start < len(rows):
        stop = max(int(np.searchsorted(cost, (cost[start - 1] if start else 0) + 1.5e8, side="right")), start + 1)
        R = rows[start:stop]
        for M, out in ((A, c), (Aw, ra)):
            Sb = (A[R] @ M.T).tocsr()
            Sb.sort_indices()
            sk = np.repeat(R.astype(np.int64), np.diff(Sb.indptr)) * n + Sb.indices
            lo, hi = np.searchsorted(keys, [int(R[0]) * n, (int(R[-1]) + 1) * n])
            q = keys[lo:hi]
            if len(sk):
                pos = np.minimum(np.searchsorted(sk, q), len(sk) - 1)
                out[lo:hi] = np.where(sk[pos] == q, Sb.data[pos], 0)
        start = stop
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    c, ra = c[inv], ra[inv]
    adj = S.pair_is_edge
    I = c + 2 * adj
    su, sv = d[pu] + 1.0, d[pv] + 1.0
    return {"I": I, "U": su + sv - I, "su": su, "sv": sv, "ra": ra}


def task_hc(net, setting, out):
    from hc_scale import ScaleScores
    e, ids = load(net, setting)
    n = len(ids)
    A = adjacency(e, n)
    stores = {m: Store(out, net, setting, m) for m in ("hc", "hc_tri", "hc_2test", "hc_econn", "linkcomm")}
    fstores = {f: Store(out, net, setting, "hcf_" + f) for f in FUNCS}
    if all(s.done() for s in list(stores.values()) + list(fstores.values())):
        return
    t0 = time.time()
    S = ScaleScores(edges=e, n=n)
    S.nodes = list(range(n))
    build = round(time.time() - t0, 1)
    for m, st in stores.items():
        if st.done():
            continue
        for t in T_GRID:
            t1 = time.time()
            cov = S.linkcomm(t) if m == "linkcomm" else S.communities(t, variant=m)
            st.add(t, *to_arrays(cov), A, ids, time.time() - t1, build_seconds=build, n_wedges=S.n_wedges)
        st.save()
    if all(s.done() for s in fstores.values()):
        return
    t0 = time.time()
    P = pair_quantities(S)
    pq = round(time.time() - t0, 1)
    jac = FUNCTIONS["jaccard"](P)[S.pair_of].min(axis=1)
    assert np.allclose(jac, S.s3), "family Jaccard differs from ScaleScores"
    saved = S.s3
    for f, st in fstores.items():
        if st.done():
            continue
        s = FUNCTIONS[f](P)[S.pair_of].min(axis=1)
        for q in Q50:
            thr = np.quantile(s, q) if q > 0 else -np.inf
            S.s3 = s
            t1 = time.time()
            cov = S.communities(thr, variant="hc")
            S.s3 = saved
            st.add(q, *to_arrays(cov), A, ids, time.time() - t1, build_seconds=build + pq,
                   n_wedges=S.n_wedges, threshold=float(thr))
        st.save()


def task_leiden(net, setting, method, gamma, out):
    import igraph as ig
    import leidenalg as la
    st = Store(out, net, setting, method, f"_{gamma}")
    if st.done():
        return
    e, ids = load(net, setting)
    n = len(ids)
    A = adjacency(e, n)
    g = ig.Graph(n=n, edges=e.tolist())
    P = la.RBConfigurationVertexPartition if method == "leiden" else la.CPMVertexPartition
    t1 = time.time()
    part = la.find_partition(g, P, resolution_parameter=gamma, n_iterations=-1, seed=1)
    st.add(gamma, *to_arrays([set(c) for c in part if len(c)]), A, ids, time.time() - t1)
    st.save()


def task_infomap(net, setting, mt, out):
    from infomap import Infomap
    st = Store(out, net, setting, "infomap", f"_{mt}")
    if st.done():
        return
    e, ids = load(net, setting)
    n = len(ids)
    A = adjacency(e, n)
    t1 = time.time()
    im = Infomap(f"--two-level --silent --seed 1 --markov-time {mt}")
    for a, b in e:
        im.add_link(int(a), int(b))
    im.run()
    groups = {}
    for node in im.tree:
        if node.is_leaf:
            groups.setdefault(node.module_id, set()).add(node.node_id)
    st.add(mt, *to_arrays(list(groups.values())), A, ids, time.time() - t1)
    st.save()


def task_simple(net, setting, out):
    """Connected components and label propagation; on induced graphs CPM and SLPA."""
    import networkx as nx
    e, ids = load(net, setting)
    n = len(ids)
    A = adjacency(e, n)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(map(tuple, e.tolist()))
    st = Store(out, net, setting, "components")
    if not st.done():
        t1 = time.time()
        st.add("", *to_arrays([set(c) for c in nx.connected_components(G)]), A, ids, time.time() - t1)
        st.save()
    st = Store(out, net, setting, "lpa")
    if not st.done():
        t1 = time.time()
        st.add("", *to_arrays([set(c) for c in nx.community.label_propagation_communities(G)]), A, ids,
               time.time() - t1)
        st.save()
    if setting != "induced":
        return
    st = Store(out, net, setting, "cpm")
    if not st.done():
        for k in (3, 4, 5, 6, 7, 8):
            t1 = time.time()
            st.add(k, *to_arrays([set(c) for c in nx.community.k_clique_communities(G, k)]), A, ids,
                   time.time() - t1)
        st.save()
    st = Store(out, net, setting, "slpa")
    if not st.done():
        from slpa import slpa_multi
        np.random.seed(1)
        t1 = time.time()
        for r, cov in slpa_multi(G, 100, [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45]):
            st.add(r, *to_arrays(cov), A, ids, time.time() - t1)
            t1 = time.time()
        st.save()


def run_task(args):
    name, a = args[0], args[1:]
    t0 = time.time()
    try:
        {"hc": task_hc, "leiden": task_leiden, "infomap": task_infomap, "simple": task_simple}[name](*a)
        return f"{time.strftime('%H:%M:%S')} done  {name} {a[:-1]} ({(time.time() - t0) / 60:.1f} min)"
    except Exception:
        with open(os.path.join(a[-1], "errors.txt"), "a") as fh:
            fh.write(f"{name} {a[:-1]}\n{traceback.format_exc()}\n")
        return f"{time.strftime('%H:%M:%S')} ERROR {name} {a[:-1]} (see errors.txt)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--only", default="amazon,dblp,youtube")
    ap.add_argument("--settings", default="full,induced")
    ap.add_argument("--out", default=os.path.join(HERE, "results_big"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tasks = []
    for net in a.only.split(","):
        for setting in a.settings.split(","):
            if not (net == "youtube" and setting == "full"):
                tasks.append(("hc", net, setting, a.out))
            tasks += [("leiden", net, setting, m, g, a.out) for m, grid in (("leiden", LEIDEN_RB),
                                                                          ("leiden_cpm", LEIDEN_CPM)) for g in grid]
            tasks += [("infomap", net, setting, mt, a.out) for mt in MARKOV]
            tasks.append(("simple", net, setting, a.out))
    # the largest graphs first, so that the long tasks start early
    tasks.sort(key=lambda t: (t[1] != "youtube", t[2] != "full", t[0] != "hc"))
    log = open(os.path.join(a.out, "log.txt"), "a")
    msg = f"{time.strftime('%H:%M:%S')} start: {len(tasks)} tasks, {a.workers} workers"
    print(msg, flush=True); log.write(msg + "\n"); log.flush()
    with Pool(a.workers, maxtasksperchild=1) as pool:
        for i, m in enumerate(pool.imap_unordered(run_task, tasks), 1):
            msg = f"[{i}/{len(tasks)}] {m}"
            print(msg, flush=True); log.write(msg + "\n"); log.flush()
    msg = f"{time.strftime('%H:%M:%S')} ALL DONE. Please send back the whole results_big folder."
    print(msg, flush=True); log.write(msg + "\n")


if __name__ == "__main__":
    main()
