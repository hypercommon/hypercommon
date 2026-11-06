"""
Community detection on SNAP Amazon and DBLP, main runs.

Changes from real_detection_v3.py: the method set -- HyperCommon and its
ablations (triangles only, two-test, edge-connected), link communities cut at every
threshold, Leiden with RB modularity and with the Constant Potts Model, Infomap over
Markov times, label propagation, and on the induced subnetworks clique percolation and
SLPA -- and every cover scored with `eq`, its overlapping modularity on the graph
detected on, so that every method's label-free choice (largest eq) is read from the
same rows. Seconds per cover and the HyperCommon build time and wedge count are
recorded. Scoring as in v3: detected communities restricted to the nodes of the
top-5000 communities; every nonempty restricted community kept (f1, coverage, ...) or
only those of at least three nodes (*_ge3); ONMI on each method's best-F1 row and
largest-eq row.

    python real_detection_v4.py <dblp|amazon|youtube> <full|induced>
Rows go to results/real_detection_v4.csv (appended; reruns skip methods already done).
"""
from __future__ import annotations

import csv
import os
import sys
import time
import warnings

import networkx as nx
import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cover_metrics as M                                     # noqa: E402
from hc_scale import ScaleScores                              # noqa: E402
from real_detection_v3 import load, onmi_large, restrict     # noqa: E402
from hc_big import overlapping_modularity as om_arrays       # noqa: E402

OUT = os.path.join(HERE, "results", "real_detection_v4.csv")
FIELDS = ["network", "setting", "method", "param", "eq", "f1", "f1_dedup", "onmi", "coverage", "n_comms",
          "f1_ge3", "f1_dedup_ge3", "coverage_ge3", "n_comms_ge3", "seconds", "build_seconds", "n_wedges",
          "n_nodes", "n_edges"]
# YouTube full graph: HyperCommon covers computed on a 128 GB machine by youtube_job.py
YT_FULL = os.path.join(HERE, "results", "youtube_full_dl", "results")
T_GRID = [round(0.02 * i, 2) for i in range(1, 31)]
LEIDEN_RB = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
             256.0, 512.0, 1024.0, 2048.0, 4096.0]
LEIDEN_CPM = [float(f"{x:.3g}") for x in np.geomspace(1e-5, 0.5, 20)]
MARKOV = [float(f"{x:.4g}") for x in np.geomspace(0.25, 16, 13)]


def overlapping_modularity(A, idx, cover):
    """self_tune.overlapping_modularity, vectorised (hc_big; test_hc_big checks equality)."""
    node = np.fromiter((idx[u] for c in cover for u in c), np.int64)
    comm = np.repeat(np.arange(len(cover)), [len(c) for c in cover])
    return om_arrays(sp.csr_matrix(A), node, comm)


def youtube_full_covers(method):
    """(t, cover, eq, seconds) for each t of the grid, read from youtube_job.py's output."""
    import json
    meta = json.load(open(os.path.join(YT_FULL, "meta.json")))
    rows = {}
    with open(os.path.join(YT_FULL, "summary.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            if r["variant"] == method:
                rows[float(r["t"])] = r
    z = np.load(os.path.join(YT_FULL, f"covers_{method}.npz"))
    for t in T_GRID:
        node, size = z[f"t{t:.2f}_node"], z[f"t{t:.2f}_size"]
        cover = [set(b.tolist()) for b in np.split(node, np.cumsum(size)[:-1])] if len(size) else []
        yield t, cover, float(rows[t]["eq"]), float(rows[t]["seconds"]), meta


def run(net, setting):
    G, truth, truth_dedup = load(net)
    universe = set().union(*truth)
    H = G.subgraph(universe).copy() if setting == "induced" else G
    nodes = list(H.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    e = np.array([(idx[a], idx[b]) for a, b in H.edges()], dtype=np.int64)
    import scipy.sparse as sp
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                      shape=(len(nodes), len(nodes)))
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            done = {(r["network"], r["setting"], r["method"]) for r in csv.DictReader(fh)}
    new = not os.path.exists(OUT)
    fh = open(OUT, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    buf = []

    def emit(method, param, cover, secs, build="", W="", eq=None):
        cov = restrict(cover, universe, 1)
        cov3 = [c for c in cov if len(c) >= 3]
        buf.append(dict(network=net, setting=setting, method=method, param=param,
                        eq=overlapping_modularity(A, idx, cover) if eq is None else eq, seconds=round(secs, 2),
                        build_seconds=build, n_wedges=W, n_nodes=H.number_of_nodes(), n_edges=H.number_of_edges(),
                        f1=M.f1_best_match(cov, truth, universe),
                        f1_dedup=M.f1_best_match(cov, truth_dedup, universe), onmi="",
                        coverage=M.coverage(cov, truth), n_comms=len(cov),
                        f1_ge3=M.f1_best_match(cov3, truth, universe),
                        f1_dedup_ge3=M.f1_best_match(cov3, truth_dedup, universe),
                        coverage_ge3=M.coverage(cov3, truth), n_comms_ge3=len(cov3), _cover=cov))
        r = buf[-1]
        print(f"  {net} {setting} {method:10s} {str(param):8s} F1 {r['f1']:.3f} eq {r['eq']:.3f} "
              f"cov {r['coverage']:.2f} k={r['n_comms']} {secs:.1f}s", flush=True)

    def flush():
        for col in ("f1", "eq"):
            best = {}
            for r in buf:
                if r["method"] not in best or r[col] > best[r["method"]][col]:
                    best[r["method"]] = r
            for r in best.values():
                if r["onmi"] == "":
                    r["onmi"] = onmi_large(r["_cover"], truth, universe)
        for r in buf:
            w.writerow({k: v for k, v in r.items() if k != "_cover"})
        fh.flush()
        buf.clear()

    print(f"{net} {setting}: n={H.number_of_nodes()} m={H.number_of_edges()} truth={len(truth)}", flush=True)
    if (net, setting, "components") not in done:
        t1 = time.time()
        emit("components", "", [set(c) for c in nx.connected_components(H)], time.time() - t1)
        flush()

    hc_methods = ["hc", "hc_tri", "hc_2test", "hc_econn", "linkcomm"]
    if net == "youtube" and setting == "full":
        # too large for ScaleScores here; link communities admit no edge pruning (hc_big)
        for m in hc_methods[:4]:
            if (net, setting, m) in done:
                continue
            for t, cov, eq, secs, meta in youtube_full_covers(m):
                emit(m, t, cov, secs, meta["build_seconds"], meta["wedges_kept"], eq=eq)
            flush()
        hc_methods = []
    if hc_methods and any((net, setting, m) not in done for m in hc_methods):
        t0 = time.time()
        S = ScaleScores(edges=e, n=len(nodes))
        build = round(time.time() - t0, 2)
        S.nodes = nodes
        for m in hc_methods:
            if (net, setting, m) in done:
                continue
            for t in T_GRID:
                t1 = time.time()
                cov = S.linkcomm(t) if m == "linkcomm" else S.communities(t, variant=m)
                emit(m, t, cov, time.time() - t1, build, S.n_wedges)
            flush()
        del S

    if (net, setting, "leiden") not in done or (net, setting, "leiden_cpm") not in done:
        import igraph as ig
        import leidenalg as la
        g = ig.Graph(n=len(nodes), edges=e.tolist())
        for m, P, grid in (("leiden", la.RBConfigurationVertexPartition, LEIDEN_RB),
                           ("leiden_cpm", la.CPMVertexPartition, LEIDEN_CPM)):
            if (net, setting, m) in done:
                continue
            for gamma in grid:
                t1 = time.time()
                part = la.find_partition(g, P, resolution_parameter=gamma, n_iterations=-1, seed=1)
                emit(m, gamma, [{nodes[i] for i in c} for c in part if len(c)], time.time() - t1)
            flush()

    if (net, setting, "infomap") not in done:
        from infomap import Infomap
        for mt in MARKOV:
            t1 = time.time()
            im = Infomap(f"--two-level --silent --seed 1 --markov-time {mt}")
            for a, b in e:
                im.add_link(int(a), int(b))
            im.run()
            groups = {}
            for node in im.tree:
                if node.is_leaf:
                    groups.setdefault(node.module_id, set()).add(nodes[node.node_id])
            emit("infomap", mt, list(groups.values()), time.time() - t1)
        flush()

    if (net, setting, "lpa") not in done:
        t1 = time.time()
        emit("lpa", "", [set(c) for c in nx.community.label_propagation_communities(H)], time.time() - t1)
        flush()

    if setting == "induced" and (net, setting, "cpm") not in done:
        for k in (3, 4, 5, 6, 7, 8):
            t1 = time.time()
            emit("cpm", k, [set(c) for c in nx.community.k_clique_communities(H, k)], time.time() - t1)
        flush()

    if setting == "induced" and (net, setting, "slpa") not in done:
        from methods2 import _slpa_multi
        np.random.seed(1)
        t1 = time.time()
        for r, cov in _slpa_multi(H, 100, [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45]):
            emit("slpa", r, cov, time.time() - t1)
            t1 = time.time()
        flush()
    fh.close()


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
