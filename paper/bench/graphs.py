"""
Benchmark graphs, each fully determined by (family, config, p, seed).

rings  ring lattices of given sizes and even degrees (generators.ring_lattice),
       floor(o*n) inter-ring node merges (generators.overlap), then a fraction p
       of the ORIGINAL edges rewired: each keeps one endpoint (chosen at random)
       and moves the other to a uniformly random non-neighbour. Every (p, seed)
       is an independent graph -- no shared trajectory between p values.
lfr    overlapping LFR (Lancichinetti & Fortunato 2009), loaded from files made
       by the reference generator (see lfr.py).
"""

from __future__ import annotations

import math
import os
import random
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.append(ROOT)   # after bench/, so repo-root "experiments" never shadows ours

from generators.ring_lattice import ring_lattice          # noqa: E402
from generators.overlap import apply_overlap, ground_truth_with_overlap  # noqa: E402


def _rewire_inside(G, sizes, q, rng):
    """Rewire a fraction q of each ring's edges INSIDE the ring: each keeps one endpoint
    and moves the other to a uniformly random non-neighbour of the same ring. q = 0 keeps
    the lattice; q = 1 leaves each ring a random graph with the same number of edges
    (short diameter, clustering falling with q)."""
    off = 0
    for s in sizes:
        ring = list(range(off, off + s))
        edges = [(u, v) for u in ring for v in G.adj[u] if off <= v < off + s and u < v]
        rng.shuffle(edges)
        for u, v in edges[:math.floor(q * len(edges) + 1e-9)]:
            keep, move = (u, v) if rng.random() < 0.5 else (v, u)
            cand = [w for w in ring if w != keep and w not in G.adj[keep]]
            if not cand:
                continue
            G.remove_edge(keep, move)
            G.add_edge(keep, rng.choice(cand))
        off += s


def _isolated_overlap(G, sizes, ks, overlap, rng):
    """Like generators.overlap.apply_overlap, but every merge is isolated in the sense of
    Theorem 2: merged positions on any one ring are more than 2k apart (cyclically), so
    no two merges share a wedge. Returns {absorbed: surviving}."""
    n = sum(sizes)
    target = math.floor(overlap * n)
    offs = [sum(sizes[:r]) for r in range(len(sizes))]
    ring_of = [r for r, s in enumerate(sizes) for _ in range(s)]
    taken = {r: [] for r in range(len(sizes))}          # merged positions per ring

    def free(x):
        r = ring_of[x]
        pos, s, k = x - offs[r], sizes[r], ks[r]
        return all(min(abs(pos - q), s - abs(pos - q)) > 2 * k for q in taken[r])

    merged = {}
    for _ in range(target * 200):
        if len(merged) >= target:
            break
        u, v = rng.randrange(n), rng.randrange(n)
        if ring_of[u] == ring_of[v] or not (free(u) and free(v)):
            continue
        for w in list(G.neighbors(v)):
            if w != u:
                G.add_edge(u, w)
        G.remove_node(v)
        merged[v] = u
        taken[ring_of[u]].append(u - offs[ring_of[u]])
        taken[ring_of[v]].append(v - offs[ring_of[v]])
    return merged


def _swap_inside(G, sizes, q, rng):
    """Degree-preserving counterpart of _rewire_inside: double-edge swaps inside each
    ring ({a,b},{c,d} -> {a,d},{c,b}, rejected when they would create a loop or a
    multi-edge) until a fraction q of the ring's edges has been moved. Every node keeps
    its degree; q = 1 leaves each ring close to a random regular graph."""
    off = 0
    for s in sizes:
        ring = range(off, off + s)
        edges = [(u, v) for u in ring for v in G.adj[u] if off <= v < off + s and u < v]
        target, moved, tries = math.floor(q * len(edges) / 2 + 1e-9), 0, 0
        while moved < target and tries < 100 * max(target, 1):
            tries += 1
            (a, b), (c, d) = rng.sample(edges, 2)
            if rng.random() < 0.5:
                c, d = d, c
            if len({a, b, c, d}) < 4 or G.has_edge(a, d) or G.has_edge(c, b):
                continue
            G.remove_edge(a, b)
            G.remove_edge(c, d)
            G.add_edge(a, d)
            G.add_edge(c, b)
            edges.remove((a, b) if a < b else (b, a))
            edges.remove((c, d) if c < d else (d, c))
            edges += [tuple(sorted((a, d))), tuple(sorted((c, b)))]
            moved += 1
        off += s


def _hub_lattice(sizes, ks, gamma, kmax_factor, rng):
    """Ring lattices with hubs: node i of ring r draws its own even degree k_i from a
    discrete power law P(k) ~ k^-gamma on {k_r, k_r + 2, ..., kmax}, kmax =
    min(kmax_factor * k_r, s_r - 2), and is joined to its k_i/2 nearest neighbours on
    either side; an edge exists when either endpoint reaches the other, so every node
    keeps its lattice neighbourhood and hubs also reach further along the ring."""
    G = nx.Graph()
    G.add_nodes_from(range(sum(sizes)))
    off = 0
    for s, k in zip(sizes, ks):
        support = list(range(k, min(kmax_factor * k, s - 2) + 1, 2))
        w = [x ** -gamma for x in support]
        for i in range(s):
            half = rng.choices(support, weights=w)[0] // 2
            for step in range(1, half + 1):
                G.add_edge(off + i, off + (i + step) % s)
                G.add_edge(off + i, off + (i - step) % s)
        off += s
    return G


def _rgg_communities(sizes, ks, rng):
    """Each community a random geometric graph on the unit torus with mean degree k:
    s points, joined when closer than r with pi r^2 (s - 1) = k. Clustered (about 0.59
    in two dimensions) but not a lattice."""
    G = nx.Graph()
    G.add_nodes_from(range(sum(sizes)))
    off = 0
    for s, k in zip(sizes, ks):
        r2 = k / (math.pi * (s - 1))
        pts = [(rng.random(), rng.random()) for _ in range(s)]
        for i in range(s):
            for j in range(i + 1, s):
                dx = abs(pts[i][0] - pts[j][0])
                dy = abs(pts[i][1] - pts[j][1])
                dx, dy = min(dx, 1 - dx), min(dy, 1 - dy)
                if dx * dx + dy * dy < r2:
                    G.add_edge(off + i, off + j)
        off += s
    return G


def _split_overlap(G, sizes, overlap, rng):
    """Like generators.overlap.apply_overlap, but the merged node SPLITS its degree:
    before the merge each of the two nodes drops a random half of its edges, so the
    survivor keeps about k/2 edges in each ring instead of k in each. The case least
    favourable to commonality, as for LFR's overlapping nodes."""
    from generators.ring_lattice import ring_of_node
    n = sum(sizes)
    target = math.floor(overlap * n)
    node_ring = ring_of_node(sizes)
    merged, used = {}, set()
    for _ in range(target * 20):
        if len(merged) >= target:
            break
        u, v = rng.randrange(n), rng.randrange(n)
        if u == v or node_ring[u] == node_ring[v] or u in used or v in used:
            continue
        if not (G.has_node(u) and G.has_node(v)):
            continue
        for x in (u, v):
            nb = sorted(G.neighbors(x))
            for y in rng.sample(nb, len(nb) // 2):
                G.remove_edge(x, y)
        for w in list(G.neighbors(v)):
            if w != u:
                G.add_edge(u, w)
        G.remove_node(v)
        used |= {u, v}
        merged[v] = u
    return merged


def rings(sizes, ks, overlap, p, seed, rewire="uniform", internal=0.0, isolated=False,
          internal_mode="endpoint", split=False, hubs=None, geometry="lattice"):
    """
    rewire="uniform"       the moved endpoint goes to a uniformly random non-neighbour
    rewire="preferential"  ... to node w with probability proportional to deg(w)+1,
                           using CURRENT degrees (the degree-biased randomisation of
                           the thesis); rich-get-richer, but mild from a lattice start
    rewire="preferential2" ... proportional to (deg(w)+1)^2: superlinear, grows stars
    internal = q           first rewire a fraction q of each ring's edges inside the ring
                           (community geometry: lattice at q = 0, random graph at q = 1)
    internal_mode="swap"   ... by degree-preserving double-edge swaps instead of moving
                           one endpoint (every node keeps its degree)
    isolated = True        merges placed more than 2k apart on every ring (Theorem 2)
    split = True           merged nodes split their degree between the two rings
    hubs = (gamma, f)      ring lattices with power-law node degrees (_hub_lattice)
    geometry="rgg"         communities are random geometric graphs, not rings
    """
    tag = "" if rewire == "uniform" else f"|{rewire}"
    tag += f"|in{internal}" if internal else ""
    tag += "|iso" if isolated else ""
    tag += "|swap" if internal_mode == "swap" else ""
    tag += "|split" if split else ""
    tag += f"|hubs{hubs}" if hubs else ""
    tag += "|rgg" if geometry == "rgg" else ""
    rng = random.Random(f"rings|{sizes}|{ks}|{overlap}|{p}|{seed}{tag}")
    if geometry == "rgg":
        G = _rgg_communities(list(sizes), list(ks), rng)
    elif hubs:
        G = _hub_lattice(list(sizes), list(ks), hubs[0], hubs[1], rng)
    else:
        G = ring_lattice(list(sizes), list(ks))
    if internal:
        (_swap_inside if internal_mode == "swap" else _rewire_inside)(G, list(sizes), internal, rng)
    if isolated:
        merged = _isolated_overlap(G, list(sizes), list(ks), overlap, rng)
    elif split:
        merged = _split_overlap(G, list(sizes), overlap, rng)
    else:
        merged = apply_overlap(G, list(sizes), overlap, rng)
    truth = ground_truth_with_overlap(list(sizes), merged)
    edges = list(G.edges())
    rng.shuffle(edges)
    k = math.floor(p * len(edges) + 1e-9)
    nodes = list(G.nodes())
    alpha = {"uniform": 0.0, "preferential": 1.0, "preferential2": 2.0}[rewire]
    weight = [(G.degree(u) + 1) ** alpha for u in nodes] if alpha else None
    pos = {u: i for i, u in enumerate(nodes)}
    for u, v in edges[:k]:
        keep, move = (u, v) if rng.random() < 0.5 else (v, u)
        nbrs = G.adj[keep]
        while True:                       # rejection sample a non-neighbour
            w = rng.choice(nodes) if weight is None else rng.choices(nodes, weights=weight)[0]
            if w != keep and w not in nbrs:
                break
        G.remove_edge(keep, move)
        G.add_edge(keep, w)
        if weight is not None:            # current degrees: rich get richer
            for x in (move, w):
                weight[pos[x]] = (G.degree(x) + 1) ** alpha
    return G, [set(c) for c in truth]


# --------------------------------------------------------------------------- LFR
# Reference generator: github.com/andrealancichinetti/LFRbenchmarks @ ec9a860
# (unweighted_undirected, compiled unmodified with g++ -O3 -std=c++98). It reads
# its seed from time_seed.dat in the working directory, so each graph is made in
# its own cache folder with the seed written first -- fully reproducible.

import hashlib
import shutil
import subprocess
import networkx as nx

LFR_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin", "lfr_benchmark.exe")
LFR_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lfr_cache")


def lfr(params, mu, seed):
    """params: dict with N, k, maxk, t1, t2, minc, maxc, on_frac, om."""
    N = params["N"]
    on = int(round(params["on_frac"] * N))
    args = ["-N", N, "-k", params["k"], "-maxk", params["maxk"], "-mu", mu,
            "-t1", params["t1"], "-t2", params["t2"], "-minc", params["minc"],
            "-maxc", params["maxc"], "-on", on, "-om", params["om"] if on else 0]
    if params.get("C") is not None:          # target average clustering (LFR's own rewiring)
        args += ["-C", params["C"]]
    args = [str(a) for a in args]
    key = hashlib.sha1((" ".join(args) + f"|{seed}").encode()).hexdigest()[:16]
    d = os.path.join(LFR_CACHE, key)
    net, com = os.path.join(d, "network.dat"), os.path.join(d, "community.dat")
    if not (os.path.exists(net) and os.path.exists(com)):
        tmp = d + f".tmp{os.getpid()}"
        os.makedirs(tmp, exist_ok=True)
        with open(os.path.join(tmp, "time_seed.dat"), "w") as fh:
            fh.write(str(seed + 1))
        subprocess.run([LFR_BIN] + args, cwd=tmp, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True, timeout=600)
        with open(os.path.join(tmp, "args.txt"), "w") as fh:
            fh.write(" ".join(args) + f"\nseed {seed + 1}\n")
        try:
            os.replace(tmp, d)
        except OSError:              # another worker won the race
            shutil.rmtree(tmp, ignore_errors=True)
    G = nx.Graph()
    with open(net) as fh:
        for line in fh:
            a, b = line.split()[:2]
            if a != b:
                G.add_edge(int(a), int(b))
    comms = {}
    with open(com) as fh:
        for line in fh:
            parts = line.split()
            u = int(parts[0])
            G.add_node(u)
            for c in parts[1:]:
                comms.setdefault(int(c), set()).add(u)
    return G, list(comms.values())


# ------------------------------------------------------------------------- nulls
def null_graph(kind, seed):
    """Clustered graphs WITHOUT planted communities (truth: one community of all
    nodes): 'rgg', one random geometric graph on the unit torus (n = 2000, mean degree
    8); 'holme_kim', the Holme-Kim power-law graph with triadic closure (n = 2000,
    m = 4, triangle probability 0.9)."""
    rng = random.Random(f"null|{kind}|{seed}")
    if kind == "rgg":
        G = _rgg_communities([2000], [8], rng)
    elif kind == "holme_kim":
        G = nx.powerlaw_cluster_graph(2000, 4, 0.9, seed=rng.randrange(2**31))
    else:
        raise ValueError(kind)
    return G, [set(G.nodes())]
