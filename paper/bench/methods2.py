"""
Method set and grids of the second benchmark round (run2.py), fixed before any of
its results were seen. Changes from methods.py:

  leiden      Reichardt-Bornholdt modularity on 30 log-spaced resolutions
              (0.0625 .. 32) instead of 6, best of 3 restarts per resolution,
              iterated to convergence
  leiden_cpm  Leiden with the Constant Potts Model (resolution-limit-free;
              Traag, Van Dooren & Nesterov 2011), 30 log-spaced resolutions
              (0.002 .. 0.5), same restarts
  infomap     two-level Infomap over 13 Markov times (0.25 .. 16, including 1)
  linkcomm    link communities (Ahn, Bagrow & Lehmann 2010) cut at every value of
              the HyperCommon threshold grid instead of only at maximum partition
              density; the link similarity is the same closed-neighbourhood Jaccard
  slpa        100 iterations everywhere, r in {0.01 .. 0.45}
  hc_econn    ablation: HyperCommon's triples, joined only through shared EDGES
              (not through shared non-adjacent pairs)
  cpm, demon, angel, scan, scan_j: their extended grids from the first round.

Every method's every cover is scored, so the oracle, frozen and label-free
protocols are all read from the same rows (analysis2.py).
"""

from __future__ import annotations

import random
import warnings

import numpy as np

from methods import T_GRID, _scan, _sets, _structural_similarity

warnings.filterwarnings("ignore")


def _geom(lo, hi, n):
    return [float(f"{x:.4g}") for x in np.geomspace(lo, hi, n)]


GRIDS2 = {
    "hc": T_GRID, "hc_2test": T_GRID, "hc_tri": T_GRID, "hc_econn": T_GRID,
    "linkcomm": T_GRID,
    "cpm": [3, 4, 5, 6, 7],
    "leiden": _geom(1 / 16, 32, 30),
    "leiden_cpm": _geom(0.002, 0.5, 30),
    "infomap": _geom(0.25, 16, 13),
    "lpa": [None],
    "walktrap": [None],
    "slpa": [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45],
    "demon": [0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
    "angel": [0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
    "scan": [(e, m) for e in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8) for m in (2, 3, 5)],
    "scan_j": [(e, m) for e in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6) for m in (2, 3, 5)],
}
METHODS2 = ["hc", "hc_2test", "hc_tri", "hc_econn", "linkcomm", "cpm", "leiden", "leiden_cpm",
            "infomap", "lpa", "walktrap", "slpa", "demon", "angel", "scan", "scan_j"]
HC_LIKE = ("hc", "hc_2test", "hc_tri", "hc_econn", "linkcomm", "linkcomm_pd")
# Native defaults (each method's own choice, no sweep), for the protocol that does not
# select by overlapping modularity: Leiden-RB at gamma = 1; Leiden-CPM at gamma equal
# to the edge density 2m/(n(n-1)) (the Erdos-Renyi null of the CPM); link communities
# at the cut of largest partition density over the threshold grid (Ahn et al.'s rule).
NATIVE = ["leiden_g1", "leiden_cpm_dens", "linkcomm_pd"]
RESTARTS = 3


def _leiden(G, kind, gammas, seed):
    import igraph as ig
    import leidenalg as la
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    g = ig.Graph(n=len(nodes), edges=[(idx[u], idx[v]) for u, v in G.edges()])
    P = la.RBConfigurationVertexPartition if kind == "rb" else la.CPMVertexPartition
    for gamma in gammas:
        best = None
        for r in range(RESTARTS):
            part = la.find_partition(g, P, resolution_parameter=gamma, n_iterations=-1,
                                     seed=(seed * RESTARTS + r) % (2**31))
            if best is None or part.quality() > best.quality():
                best = part
        yield gamma, [{nodes[i] for i in c} for c in best if len(c) > 0]


def _infomap(G, markov_times, seed):
    from infomap import Infomap
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    for mt in markov_times:
        im = Infomap(f"--two-level --silent --seed {seed % (2**31) + 1} --markov-time {mt}")
        for u, v in G.edges():
            im.add_link(idx[u], idx[v])
        im.run()
        groups = {}
        for node in im.tree:
            if node.is_leaf:
                groups.setdefault(node.module_id, set()).add(nodes[node.node_id])
        yield mt, [c for c in groups.values() if len(c) > 0]


def _slpa_multi(G, T, rs):
    """SLPA (Xie, Szymanski & Liu 2011) as in cdlib's slpa_nx, with the label-propagation
    stage run once and its post-processing threshold r applied for every r of the grid
    (the threshold only filters the final memories)."""
    from collections import defaultdict
    memory = {i: {i: 1} for i in G.nodes()}
    for _ in range(T):
        order = list(G.nodes())
        np.random.shuffle(order)
        for listener in order:
            speakers = G[listener].keys()
            if len(speakers) == 0:
                continue
            labels = defaultdict(int)
            for speaker in speakers:
                mem = memory[speaker]
                total = float(sum(mem.values()))
                labels[list(mem.keys())[np.random.multinomial(1, [f / total for f in mem.values()]).argmax()]] += 1
            acc = max(labels, key=labels.get)
            memory[listener][acc] = memory[listener].get(acc, 0) + 1
    for r in rs:
        comms = {}
        for node, mem in memory.items():
            for label, freq in mem.items():
                if freq / float(T + 1) >= r:
                    comms.setdefault(label, set()).add(node)
        keys = list(comms)
        nested = set()
        for i, a in enumerate(keys[:-1]):
            for b in keys[i + 1:]:
                if comms[a] <= comms[b]:
                    nested.add(a)
                elif comms[a] >= comms[b]:
                    nested.add(b)
        yield r, [c for k, c in comms.items() if k not in nested and len(c) > 0]


def run_method2(name, G, ts, seed):
    """Yield (param, cover) over the method's whole grid. `ts`: TripleScores of G."""
    if name in ("hc", "hc_2test", "hc_tri", "hc_econn"):
        for t in GRIDS2[name]:
            yield t, ts.communities(name, t)
        return
    if name == "linkcomm":
        for s in GRIDS2[name]:
            yield s, ts.linkcomm(s)
        return
    if name == "linkcomm_pd":
        dens = [ts.linkcomm_density(s) for s in GRIDS2["linkcomm"]]
        s = GRIDS2["linkcomm"][int(np.argmax(dens))]
        yield s, ts.linkcomm(s)
        return
    if name == "leiden_g1":
        yield from _leiden(G, "rb", [1.0], seed)
        return
    if name == "leiden_cpm_dens":
        n, m = G.number_of_nodes(), G.number_of_edges()
        yield from _leiden(G, "cpm", [2 * m / (n * (n - 1))], seed)
        return
    if name in ("leiden", "leiden_cpm"):
        yield from _leiden(G, "rb" if name == "leiden" else "cpm", GRIDS2[name], seed)
        return
    if name == "infomap":
        yield from _infomap(G, GRIDS2[name], seed)
        return
    if name == "oslom":                 # the OSLOM2 binary, as run_oslom.py (not in METHODS2)
        import shutil
        import tempfile
        from run_oslom import T_OSLOM, oslom_cover
        workdir = tempfile.mkdtemp(prefix="oslom_")
        try:
            for t in T_OSLOM:
                yield t, oslom_cover(G, t, seed, workdir)[0]
        finally:
            shutil.rmtree(workdir, ignore_errors=True)
        return

    from cdlib import algorithms as A
    random.seed(seed)
    np.random.seed(seed % (2**32))
    if name == "cpm":
        # the output of cdlib's kclique (networkx k_clique_communities), in bounded memory
        from cpm_fast import k_clique_covers
        covers = k_clique_covers(G, GRIDS2[name])
        for k in GRIDS2[name]:
            yield k, [set(c) for c in covers[k]]
    elif name == "lpa":
        yield None, _sets(A.label_propagation(G))
    elif name == "walktrap":
        yield None, _sets(A.walktrap(G))
    elif name == "slpa":
        yield from _slpa_multi(G, 100, GRIDS2[name])
    elif name == "demon":
        for e in GRIDS2[name]:
            yield e, _sets(A.demon(G, epsilon=e, min_com_size=3))
    elif name == "angel":
        for th in GRIDS2[name]:
            yield th, _sets(A.angel(G, threshold=th, min_community_size=3))
    elif name in ("scan", "scan_j"):
        sim = _structural_similarity(G, cosine=(name == "scan"))
        for e, m in GRIDS2[name]:
            yield (e, m), _scan(G, sim, e, m)
    else:
        raise ValueError(name)
