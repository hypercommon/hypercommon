"""
Community detection on real networks with ground truth (SNAP DBLP, Amazon),
following Yang & Leskovec: ground truth = the top-5000 communities, evaluated on
the subnetwork induced by their nodes ("induced"). A sensitivity setting detects
on the FULL graph and evaluates the detected communities restricted to those
nodes ("full").

Every parameterised method is scored over its whole grid; we report its best
parameter by F1 (oracle, same budget for all) and its label-free choice:
HyperCommon by maximal overlapping modularity of its own cover (self_tune.py),
Leiden at the standard resolution gamma = 1. Metrics: best-match F1, ONMI
(universe = evaluated nodes), coverage, number of communities, seconds.

    python real_detection.py <dblp|amazon> <induced|full>
Rows go to results/real_detection.csv (appended; reruns skip done methods).
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
import time
import warnings

import networkx as nx
import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cover_metrics as M                       # noqa: E402
from hc_scale import ScaleScores                # noqa: E402
from self_tune import overlapping_modularity    # noqa: E402

SNAP = os.path.join(HERE, "..", "analysis", "snap")
OUT = os.path.join(HERE, "results", "real_detection.csv")
FIELDS = ["network", "setting", "method", "param", "label_free", "f1", "onmi", "coverage",
          "n_comms", "seconds", "n_nodes", "n_edges"]
T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]


def load(net):
    e = np.loadtxt(gzip.open(os.path.join(SNAP, f"com-{net}.ungraph.txt.gz"), "rt"), comments="#", dtype=np.int64)
    G = nx.Graph()
    G.add_edges_from(map(tuple, e))
    truth = []
    with gzip.open(os.path.join(SNAP, f"com-{net}.top5000.cmty.txt.gz"), "rt") as fh:
        for line in fh:
            c = {int(x) for x in line.split()}
            if len(c) >= 3:
                truth.append(c)
    return G, truth


def _cond_entropy_chunked(Ma, Mb, sa, sb, n, chunk=1500):
    """cover_metrics._cond_entropy evaluated in row blocks, so the k_a x k_b
    contingency table is never held in memory at once."""
    h = M._h
    hy = h(sb / n) + h(1 - sb / n)
    total = 0.0
    MaT = Ma.T.tocsr()
    for r0 in range(0, Ma.shape[1], chunk):
        I = (MaT[r0:r0 + chunk] @ Mb).toarray().astype(float)
        sa_c = sa[r0:r0 + chunk][:, None]
        a = (n - sa_c - sb[None, :] + I) / n
        b = (sb[None, :] - I) / n
        c = (sa_c - I) / n
        d = I / n
        ok = h(a) + h(d) > h(b) + h(c)
        hx = h(sa_c / n) + h(1 - sa_c / n)
        cond = np.where(ok, h(a) + h(b) + h(c) + h(d) - hy[None, :], hx)
        total += float(cond.min(axis=1).sum())
    return total


def onmi_large(detected, truth, nodes):
    """Same quantity as cover_metrics.onmi (McDaid, max normalisation, uncovered
    nodes as singletons), memory-bounded."""
    nodes = list(nodes)
    n = len(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    covered = set().union(*detected) if detected else set()
    det = [set(c) for c in detected] + [{u} for u in nodes if u not in covered]
    Ma, Mb = M._membership(det, index, n), M._membership(truth, index, n)
    sa = np.asarray(Ma.sum(axis=0)).ravel().astype(float)
    sb = np.asarray(Mb.sum(axis=0)).ravel().astype(float)
    hxy = _cond_entropy_chunked(Ma, Mb, sa, sb, n)
    hyx = _cond_entropy_chunked(Mb, Ma, sb, sa, n)
    hx = float((M._h(sa / n) + M._h(1 - sa / n)).sum())
    hy = float((M._h(sb / n) + M._h(1 - sb / n)).sum())
    return 0.5 * (hx - hxy + hy - hyx) / max(hx, hy)


def restrict(cover, universe):
    cov = [set(c) & universe for c in cover]
    return [c for c in cov if len(c) >= 3]


def run(net, setting):
    G, truth = load(net)
    universe = set().union(*truth)
    Gsub = G.subgraph(universe).copy()
    H = Gsub if setting == "induced" else G
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            done = {(r["network"], r["setting"], r["method"]) for r in csv.DictReader(fh)}
    new = not os.path.exists(OUT)
    fh = open(OUT, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()

    buffer = []

    def emit(method, param, cover, secs, label_free=""):
        cov = restrict(cover, universe)
        buffer.append(dict(network=net, setting=setting, method=method, param=param, label_free=label_free,
                           seconds=round(secs, 2), n_nodes=H.number_of_nodes(), n_edges=H.number_of_edges(),
                           f1=M.f1_best_match(cov, truth, universe), onmi="",
                           coverage=M.coverage(cov, truth), n_comms=len(cov), _cover=cov))
        r = buffer[-1]
        print(f"  {net} {setting} {method:9s} {str(param):6s} F1 {r['f1']:.3f} cov {r['coverage']:.2f} "
              f"k={r['n_comms']} {secs:.0f}s {label_free}", flush=True)

    def flush_method():
        """ONMI for the best-F1 row of each method and for label-free rows, then write."""
        best = {}
        for r in buffer:
            if r["method"] not in best or r["f1"] > best[r["method"]]["f1"]:
                best[r["method"]] = r
        for r in buffer:
            if r is best[r["method"]] or r["label_free"]:
                r["onmi"] = onmi_large(r["_cover"], truth, universe)
                print(f"    ONMI {r['method']} {r['param']} {r['label_free'] or 'best-F1'}: {r['onmi']:.3f}", flush=True)
            w.writerow({k: v for k, v in r.items() if k != "_cover"})
        fh.flush()
        buffer.clear()

    print(f"{net} {setting}: n={H.number_of_nodes()} m={H.number_of_edges()} truth={len(truth)}", flush=True)

    if (net, setting, "hc") not in done:
        nodes = list(H.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        e = np.array([(idx[a], idx[b]) for a, b in H.edges()], dtype=np.int64)
        t0 = time.time()
        S = ScaleScores(edges=e, n=len(nodes))
        build = time.time() - t0
        S.nodes = nodes
        import scipy.sparse as sp
        A = S.A.astype(float)
        for variant in ("hc", "hc_tri"):
            best_q, best_t = -np.inf, None
            covers = {}
            for t in T_GRID:
                t1 = time.time()
                cov = S.communities(t, triangles_only=(variant == "hc_tri"))
                secs = time.time() - t1 + build
                covers[t] = cov
                q = overlapping_modularity(A, idx, cov)
                if q > best_q:
                    best_q, best_t = q, t
                emit(variant, t, cov, secs)
            # label-free pick, scored again under its own row
            emit(variant, best_t, covers[best_t], 0.0, label_free="max_overlapping_modularity")
            del covers
        flush_method()

    if (net, setting, "leiden") not in done:
        import igraph as ig
        import leidenalg as la
        nodes = list(H.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        g = ig.Graph(n=len(nodes), edges=[(idx[a], idx[b]) for a, b in H.edges()])
        for gamma in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
            t1 = time.time()
            part = la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter=gamma, seed=1)
            cov = [{nodes[i] for i in c} for c in part]
            emit("leiden", gamma, cov, time.time() - t1, label_free="default" if gamma == 1.0 else "")
        flush_method()

    if (net, setting, "infomap") not in done:
        from infomap import Infomap
        nodes = list(H.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        t1 = time.time()
        im = Infomap("--two-level --silent --seed 1")
        for a, b in H.edges():
            im.add_link(idx[a], idx[b])
        im.run()
        groups = {}
        for node in im.tree:
            if node.is_leaf:
                groups.setdefault(node.module_id, set()).add(nodes[node.node_id])
        emit("infomap", "", list(groups.values()), time.time() - t1, label_free="default")
        flush_method()

    if (net, setting, "lpa") not in done:
        t1 = time.time()
        cov = [set(c) for c in nx.community.label_propagation_communities(H)]
        emit("lpa", "", cov, time.time() - t1, label_free="default")
        flush_method()

    if setting == "induced" and (net, setting, "cpm") not in done:
        for k in (3, 4, 5):
            t1 = time.time()
            cov = [set(c) for c in nx.community.k_clique_communities(H, k)]
            emit("cpm", k, cov, time.time() - t1)
        flush_method()

    if setting == "induced" and (net, setting, "slpa") not in done:
        from cdlib import algorithms as A
        for r in (0.05, 0.1, 0.2, 0.3, 0.45):
            t1 = time.time()
            cov = [set(c) for c in A.slpa(H, t=21, r=r).communities]
            emit("slpa", r, cov, time.time() - t1)
        flush_method()
    fh.close()


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
