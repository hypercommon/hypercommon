"""Build sliding-window graphs from the Time Series Segmentation Benchmark (TSSB),
with the annotated segments as ground-truth groups.

Source data (read as plain text, never executed or imported):
  TSSB, Ermshaus, Schaefer, Leser, github.com/ermshaua/time-series-segmentation-benchmark
  (BSD-3-Clause), pinned to commit COMMIT below. Files used: tssb/datasets/desc.txt
  (name, window size, change points) and tssb/datasets/<name>.txt (one value per line).

Usage (run from paper/, with an isolated interpreter; downloads are untrusted data):
  python -I analysis/tssb_build.py fetch --raw D:/tmp/tssb_raw/repo
  python -I analysis/tssb_build.py build --raw D:/tmp/tssb_raw/repo --out analysis/data_tssb
  python -I analysis/tssb_build.py table --out analysis/data_tssb     # README table

`fetch` downloads desc.txt, properties.txt, LICENSE, README.md and the 75 series files at
the pinned commit from raw.githubusercontent.com into --raw (skips files already there).

Rules (fixed before any detection was run; see data_tssb/README.md):
  * Series: those with 2-9 annotated segments (69 of 75); a seeded subset of N_SERIES = 24
    of them (random.Random(SEED).sample over the sorted names), to keep <= 120 folders.
  * Nodes: sliding windows of the benchmark window size w, stride s = max(1, ceil(len/3000)),
    window i covers samples [i*s, i*s + w), i = 0..n-1, n = floor((len - w)/s) + 1.
    Each window is z-normalised (subtract its mean, divide by its population std; a window
    with std <= 1e-8 becomes the zero vector).
  * Distance: squared Euclidean distance between z-normalised windows, rounded to 9 decimals
    so that floating-point near-ties become exact ties.
  * knn graph: each window i keeps the k windows j != i with the smallest distance (ties:
    smaller time gap |i-j|, then smaller j); union over i, unweighted. No exclusion zone.
  * band graph: as knn, but candidates are only j with 0 < |i-j| <= W,
    W = max(ceil(5*w/s), 2*k) windows.
  * null control (named after the band graph): each window keeps its k nearest windows in
    time (|i-j|, ties to smaller j); union over i. Same nodes and groups.
  * Groups: segment t covers samples [c_{t-1}, c_t) with c_0 = 0, c_T = len; window i is in
    segment t if it shares at least one sample with it (i*s < c_t and i*s + w > c_{t-1}).
    Windows straddling a change point are in both segments. Groups < 3 nodes and exact
    duplicate node sets are dropped.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
import urllib.request

import numpy as np

REPO = "ermshaua/time-series-segmentation-benchmark"
COMMIT = "e0c7a1f6cee6874e60cbacd3b8601f2de259606e"
RAW_URL = "https://raw.githubusercontent.com/" + REPO + "/" + COMMIT + "/{}"
UA = "tssb-graph-builder/0.1 (python urllib)"
SEED = 2026
N_SERIES = 24
MIN_SEG, MAX_SEG = 2, 9
MAX_N = 3000
KS_KNN = (8,)
KS_BAND = (8, 16)
MIN_GROUP = 3


# ----------------------------------------------------------------------------- data
def fetch(raw):
    def get(rel):
        dst = os.path.join(raw, *rel.split("/"))
        if os.path.exists(dst):
            return dst
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        req = urllib.request.Request(RAW_URL.format(rel), headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        with open(dst, "wb") as fh:
            fh.write(data)
        time.sleep(0.5)
        return dst

    for rel in ("LICENSE", "README.md", "tssb/datasets/desc.txt", "tssb/datasets/properties.txt"):
        get(rel)
    for name, _, _ in read_desc(raw):
        get(f"tssb/datasets/{name}.txt")
    print("fetched into", raw)


def read_desc(raw):
    out = []
    with open(os.path.join(raw, "tssb", "datasets", "desc.txt"), encoding="ascii") as fh:
        for line in fh:
            p = line.strip().split(",")
            if not p[0]:
                continue
            out.append((p[0], int(p[1]), [int(c) for c in p[2:]]))
    return out


def read_series(raw, name):
    vals = []
    with open(os.path.join(raw, "tssb", "datasets", name + ".txt"), encoding="ascii") as fh:
        for line in fh:
            line = line.strip()
            if line:
                vals.append(float(line))
    x = np.asarray(vals, dtype=np.float64)
    if not np.all(np.isfinite(x)):
        raise ValueError(f"{name}: non-finite values")
    return x


def select(desc):
    elig = sorted(name for name, _, cps in desc if MIN_SEG <= len(cps) + 1 <= MAX_SEG)
    chosen = sorted(random.Random(SEED).sample(elig, N_SERIES))
    return elig, chosen


# ----------------------------------------------------------------------------- graph
def windows(x, w, s):
    n = (len(x) - w) // s + 1
    starts = np.arange(n) * s
    X = np.lib.stride_tricks.sliding_window_view(x, w)[starts].astype(np.float64)
    mu = X.mean(1, keepdims=True)
    sd = X.std(1, keepdims=True)
    flat = sd[:, 0] <= 1e-8
    Z = np.where(sd > 1e-8, (X - mu) / np.where(sd > 1e-8, sd, 1.0), 0.0)
    return starts, Z, int(flat.sum())


def sqdist(Z):
    g = Z @ Z.T
    sq = np.diag(g)
    d = sq[:, None] + sq[None, :] - 2.0 * g
    np.maximum(d, 0.0, out=d)
    return np.round(d, 9)


def topk_edges(n, k, key_fn):
    """Union over i of the first k candidates j in the order given by key_fn(i)."""
    E = set()
    for i in range(n):
        for j in key_fn(i)[:k]:
            j = int(j)
            E.add((i, j) if i < j else (j, i))
    return sorted(E)


def knn_edges(D, k, W=None):
    n = D.shape[0]
    idx = np.arange(n)

    def order(i):
        if W is None:
            cand = idx[idx != i]
        else:
            lo, hi = max(0, i - W), min(n, i + W + 1)
            cand = idx[lo:hi]
            cand = cand[cand != i]
        o = np.lexsort((cand, np.abs(cand - i), D[i, cand]))
        return cand[o]

    return topk_edges(n, k, order)


def null_edges(n, k):
    def order(i):
        lo, hi = max(0, i - k), min(n, i + k + 1)
        cand = [j for j in range(lo, hi) if j != i]
        cand.sort(key=lambda j: (abs(j - i), j))
        return cand

    return topk_edges(n, k, order)


def make_groups(starts, w, cps, length):
    bounds = [0] + list(cps) + [length]
    groups, seen, n_small, n_dup = [], set(), 0, 0
    for t in range(len(bounds) - 1):
        a, b = bounds[t], bounds[t + 1]
        mem = tuple(int(i) for i in np.nonzero((starts < b) & (starts + w > a))[0])
        if len(mem) < MIN_GROUP:
            n_small += 1
            continue
        if mem in seen:
            n_dup += 1
            continue
        seen.add(mem)
        groups.append(mem)
    return groups, n_small, n_dup


def stats(n, edges, groups, prefix=""):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    gap = np.array([b - a for a, b in edges])
    return {
        prefix + "m": len(edges),
        prefix + "mean_degree": 2 * len(edges) / n,
        prefix + "clustering": nx.average_clustering(G),
        prefix + "n_components": nx.number_connected_components(G),
        prefix + "median_gap_windows": float(np.median(gap)),
        prefix + "frac_edges_within_group": float(np.mean([bool(memb[a] & memb[b])
                                                           for a, b in edges])),
    }, gap


def random_pair_within(n, groups):
    """Fraction of all unordered node pairs that share at least one group."""
    M = np.zeros((n, len(groups)), np.float32)
    for gi, g in enumerate(groups):
        M[list(g), gi] = 1.0
    S = (M @ M.T) > 0
    return float((S.sum() - np.trace(S)) / (n * (n - 1)))


def write_graph(od, edges, groups, node_lines):
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "edges.tsv"), "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)
    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
        fh.writelines(" ".join(map(str, g)) + "\n" for g in groups)
    with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
        fh.write("id\tstart\tend\tsegments\n")
        fh.writelines(node_lines)


def cmd_build(args):
    desc = read_desc(args.raw)
    assert len(desc) == 75, len(desc)
    elig, chosen = select(desc)
    info = {name: (w, cps) for name, w, cps in desc}
    summary = []
    for name in chosen:
        w, cps = info[name]
        x = read_series(args.raw, name)
        L = len(x)
        assert all(0 < c < L for c in cps) and list(cps) == sorted(set(cps)), name
        s = max(1, math.ceil(L / MAX_N))
        starts, Z, n_flat = windows(x, w, s)
        n = len(starts)
        D = sqdist(Z)
        groups, n_small, n_dup = make_groups(starts, w, cps, L)
        sizes = [len(g) for g in groups]
        cover = np.zeros(n, int)
        bounds = [0] + list(cps) + [L]
        seg_lists = [[] for _ in range(n)]
        for t in range(len(bounds) - 1):
            for i in np.nonzero((starts < bounds[t + 1]) & (starts + w > bounds[t]))[0]:
                seg_lists[i].append(t)
        for g in groups:
            cover[list(g)] += 1
        node_lines = [f"{i}\t{int(starts[i])}\t{int(starts[i]) + w}\t"
                      f"{','.join(map(str, seg_lists[i]))}\n" for i in range(n)]
        base = dict(series=name, length=L, window=w, stride=s, n=n,
                    n_segments=len(cps) + 1, n_flat_windows=n_flat,
                    n_groups=len(groups), n_groups_dropped_small=n_small,
                    n_groups_dropped_dup=n_dup, group_size_min=min(sizes),
                    group_size_median=float(np.median(sizes)), group_size_max=max(sizes),
                    frac_nodes_in_group=float(np.mean(cover > 0)),
                    nodes_in_gt1_group=int(np.sum(cover > 1)),
                    random_pair_frac_within_group=random_pair_within(n, groups))
        nulls = {}
        for k in sorted(set(KS_KNN) | set(KS_BAND)):
            ne = null_edges(n, k)
            nulls[k] = ne
        jobs = [("knn", k, None) for k in KS_KNN] + \
               [("band", k, max(math.ceil(5 * w / s), 2 * k)) for k in KS_BAND]
        for kind, k, W in jobs:
            tag = f"{name}_{kind}_k{k}"
            edges = knn_edges(D, k, W)
            write_graph(os.path.join(args.out, tag), edges, groups, node_lines)
            row = dict(graph=tag, kind=kind, k=k, band_W=W, **base)
            st, gap = stats(n, edges, groups)
            st["frac_gap_le_k"] = float(np.mean(gap <= k))
            st["frac_windows_overlap"] = float(np.mean(gap * s < w))
            row.update(st)
            ne = nulls[k]
            nst, _ = stats(n, ne, groups, prefix="null_")
            row.update(nst)
            row["null_graph"] = f"{name}_band_k{k}_null"
            summary.append(row)
            print(json.dumps({q: row[q] for q in ("graph", "n", "m", "mean_degree",
                                                  "clustering", "n_components",
                                                  "frac_edges_within_group",
                                                  "null_frac_edges_within_group")}))
            sys.stdout.flush()
        for k in KS_BAND:
            write_graph(os.path.join(args.out, f"{name}_band_k{k}_null"), nulls[k], groups,
                        node_lines)
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(args.out, "selection.json"), "w", newline="\n") as fh:
        json.dump(dict(seed=SEED, n_series=N_SERIES, rule=f"{MIN_SEG}-{MAX_SEG} segments",
                       eligible=elig, chosen=chosen), fh, indent=1)


def cmd_table(args):
    with open(os.path.join(args.out, "summary.json"), encoding="utf-8") as fh:
        S = json.load(fh)
    print("| graph | n | w | s | W | m | mean degree | clustering | comps | groups | sizes min-med-max "
          "| within | within, control | clustering, control |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in S:
        print(f"| {r['graph']} | {r['n']} | {r['window']} | {r['stride']} | "
              f"{r['band_W'] if r['band_W'] else '-'} | {r['m']} | {r['mean_degree']:.2f} | "
              f"{r['clustering']:.3f} | {r['n_components']} | {r['n_groups']} | "
              f"{r['group_size_min']}-{r['group_size_median']:g}-{r['group_size_max']} | "
              f"{r['frac_edges_within_group']:.3f} | {r['null_frac_edges_within_group']:.3f} | "
              f"{r['null_clustering']:.3f} |")
    for kind, k in [("knn", 8), ("band", 8), ("band", 16)]:
        R = [r for r in S if r["kind"] == kind and r["k"] == k]
        f = lambda q: float(np.median([r[q] for r in R]))
        print(f"median {kind} k={k}: within {f('frac_edges_within_group'):.3f}, "
              f"control {f('null_frac_edges_within_group'):.3f}, "
              f"diff {float(np.median([r['frac_edges_within_group'] - r['null_frac_edges_within_group'] for r in R])):.4f}, "
              f"clustering {f('clustering'):.3f}, gap<=k {f('frac_gap_le_k'):.3f}, "
              f"overlap {f('frac_windows_overlap'):.3f}, max m {max(r['m'] for r in R)}, "
              f"comps>1 {sum(r['n_components'] > 1 for r in R)}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("fetch")
    a.add_argument("--raw", required=True)
    a = sub.add_parser("build")
    a.add_argument("--raw", required=True)
    a.add_argument("--out", required=True)
    a = sub.add_parser("table")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"fetch": fetch_cmd, "build": cmd_build, "table": cmd_table}[args.cmd](args)


def fetch_cmd(args):
    fetch(args.raw)


if __name__ == "__main__":
    main()
