"""Build point graphs from hand-drawn sketches (Google Quick, Draw! raw data), with the
pen strokes as ground-truth groups (long, thin chains).

Source data (downloaded by hand with curl into its own folder, never executed):
  * Quick, Draw! dataset (Google, CC BY 4.0), raw timestamped drawings,
      https://storage.googleapis.com/quickdraw_dataset/full/raw/<category>.ndjson
    Only the first 20,000,000 bytes of each category file were fetched (HTTP Range
    0-19999999); sizes and SHA-256 of these prefixes are in data_sketch/README.md.
    Files are parsed as JSON lines (json.loads), nothing else.

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/sketch_build.py build --raw <dir with <cat>.part.ndjson> --out analysis/data_sketch
  python -I analysis/sketch_build.py table --out analysis/data_sketch     # README table
(run from paper/; any working directory works, all paths are arguments).

Rules (fixed before any graph was built; no community detection was run):
  Categories: CATS (5 categories chosen because their drawings usually have many strokes).
  Scan: the first N_SCAN complete lines of each prefix file.
  Normalise: per drawing, shift x, y so the minimum is 0 and scale both by the same
    factor so that max(x range, y range) = BOX (as in the official "simplified" files).
    t (ms since the drawing's first point) is kept as is.
  Resample: per stroke, consecutive points with identical (x, y) are merged (keep the
    first); the polyline of the remaining points is resampled at arc length
    s = 0, STEP, 2*STEP, ..., floor(L/STEP)*STEP (L = stroke length; a stroke of length 0
    gives one node). x, y and t are interpolated linearly in s.
  Eligible drawing: recognized == true; >= MIN_STROKES strokes; t non-decreasing inside
    every stroke and from each stroke's last point to the next stroke's first point;
    median pen speed > 0; N_MIN <= n <= N_MAX after resampling.
  Pen speed: over all consecutive raw point pairs inside a stroke (after normalisation,
    before merging) with dt > 0, speed = |dxy| / dt (units per ms). v = median speed.
  Selection: eligible drawings of a category sorted by key_id (as integers); one
    rng = random.Random(SEED); for the categories in CATS order, rng.sample(eligible,
    PER_CAT), then sorted by key_id again.
  Node ids: in drawing order (stroke by stroke, along the stroke). x, y are rounded to
    4 decimals and t to 3 decimals; all edges are computed from the rounded values.
  Real graph <name>_k<k>: tau = v * t (time measured in "units the pen covers at its
    median speed", so 1 ms maps to v spatial units). Each node keeps its k nearest nodes
    j != i by Euclidean distance in (x, y, tau); squared distances rounded to 9 decimals,
    ties to smaller j. Undirected edge set = union over i.
  Control <name>_k<k>_null: same, in (x, y) only.
  Groups: the nodes of each stroke with >= MIN_GROUP nodes (disjoint; nodes of shorter
    strokes are in no group).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random

import numpy as np

SEED = 2026
CATS = ("bicycle", "cat", "face", "flower", "house")
PER_CAT = 4
N_SCAN = 4000
MIN_STROKES = 6
BOX = 255.0
STEP = 2.0
N_MIN, N_MAX = 300, 1500
KS = (8, 16)
MIN_GROUP = 3


# ---------------------------------------------------------------- inputs
def read_lines(path, n):
    out = []
    with open(path, "rb") as fh:
        for raw in fh:
            if not raw.endswith(b"\n"):
                break                                   # cut-off last line of the prefix
            out.append(json.loads(raw.decode("utf-8")))
            if len(out) == n:
                break
    assert len(out) == n, (path, len(out))
    return out


def parse(rec, cat):
    assert rec["word"] == cat, (rec["word"], cat)
    strokes = []
    for s in rec["drawing"]:
        assert len(s) == 3 and len(s[0]) == len(s[1]) == len(s[2]) >= 1
        strokes.append(np.array(s, dtype=np.float64).T)  # (m, 3): x, y, t
    return strokes


def normalise(strokes):
    allp = np.vstack(strokes)
    lo = allp[:, :2].min(0)
    span = (allp[:, :2].max(0) - lo).max()
    if span <= 0:
        return None
    out = []
    for s in strokes:
        s = s.copy()
        s[:, :2] = (s[:, :2] - lo) * (BOX / span)
        out.append(s)
    return out


def timing_ok(strokes):
    prev = -math.inf
    for s in strokes:
        t = s[:, 2]
        if t[0] < prev or np.any(np.diff(t) < 0):
            return False
        prev = t[-1]
    return True


def median_speed(strokes):
    sp = []
    for s in strokes:
        d = np.hypot(*np.diff(s[:, :2], axis=0).T)
        dt = np.diff(s[:, 2])
        ok = dt > 0
        sp.append(d[ok] / dt[ok])
    sp = np.concatenate(sp) if sp else np.zeros(0)
    return float(np.median(sp)) if len(sp) else 0.0


def resample(s):
    keep = np.ones(len(s), bool)
    keep[1:] = np.any(np.diff(s[:, :2], axis=0) != 0, axis=1)
    s = s[keep]
    if len(s) == 1:
        return s[:1].copy(), 0.0
    d = np.hypot(*np.diff(s[:, :2], axis=0).T)
    S = np.concatenate([[0.0], np.cumsum(d)])
    L = float(S[-1])
    q = np.arange(int(math.floor(L / STEP + 1e-9)) + 1) * STEP
    return np.column_stack([np.interp(q, S, s[:, c]) for c in range(3)]), L


def prepare(rec, cat):
    """Return None if not eligible, else a dict; 'why' says which rule failed."""
    if rec.get("recognized") is not True:
        return dict(why="unrecognized")
    strokes = parse(rec, cat)
    if len(strokes) < MIN_STROKES:
        return dict(why="few_strokes")
    strokes = normalise(strokes)
    if strokes is None or not timing_ok(strokes):
        return dict(why="bad_timing")
    v = median_speed(strokes)
    if not v > 0:
        return dict(why="zero_speed")
    pts, lens = [], []
    for s in strokes:
        p, L = resample(s)
        pts.append(p)
        lens.append(L)
    n = sum(len(p) for p in pts)
    if not N_MIN <= n <= N_MAX:
        return dict(why="n_out_of_range", n=n)
    return dict(why=None, key_id=rec["key_id"], countrycode=rec.get("countrycode"),
                timestamp=rec.get("timestamp"), raw_strokes=strokes, pts=pts,
                lens=lens, v=v, n=n)


# ---------------------------------------------------------------- edges
def knn_edges(X, k):
    n = len(X)
    E = set()
    idx = np.arange(n)
    for i in range(n):
        d2 = np.round(((X - X[i]) ** 2).sum(1), 9)
        d2[i] = np.inf
        order = np.lexsort((idx, d2))[:k]                 # last key is primary
        for j in order:
            j = int(j)
            E.add((min(i, j), max(i, j)))
    return sorted(E)


# ---------------------------------------------------------------- stats
def graph_stats(n, edges, stroke_of, pos_in_stroke, in_group):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    e = np.array(edges)
    a, b = e[:, 0], e[:, 1]
    same = stroke_of[a] == stroke_of[b]
    within = same & in_group[a] & in_group[b]
    consec = same & (np.abs(pos_in_stroke[a] - pos_in_stroke[b]) == 1)
    deg = np.bincount(e.ravel(), minlength=n)
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                degree_min=int(deg.min()), degree_max=int(deg.max()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                frac_within=float(np.mean(within)),
                frac_edges_same_stroke=float(np.mean(same)),
                frac_edges_consecutive=float(np.mean(consec)))


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


def cmd_build(args):
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(SEED)
    summary, scan = [], {}
    for cat in CATS:
        recs = read_lines(os.path.join(args.raw, f"{cat}.part.ndjson"), N_SCAN)
        preps = [prepare(r, cat) for r in recs]
        why = {}
        for p in preps:
            why[p["why"] or "eligible"] = why.get(p["why"] or "eligible", 0) + 1
        elig = sorted((p for p in preps if p["why"] is None), key=lambda p: int(p["key_id"]))
        chosen = sorted(rng.sample(elig, PER_CAT), key=lambda p: int(p["key_id"]))
        scan[cat] = dict(scanned=N_SCAN, **why,
                         chosen=[p["key_id"] for p in chosen])
        print(cat, scan[cat])
        for p in chosen:
            pts = p["pts"]
            P = np.vstack(pts)
            P[:, 0] = np.round(P[:, 0], 4)
            P[:, 1] = np.round(P[:, 1], 4)
            P[:, 2] = np.round(P[:, 2], 3)
            n = len(P)
            stroke_of = np.concatenate([np.full(len(q), si) for si, q in enumerate(pts)])
            pos = np.concatenate([np.arange(len(q)) for q in pts])
            sizes_all = np.array([len(q) for q in pts])
            groups, off = [], 0
            for q in pts:
                if len(q) >= MIN_GROUP:
                    groups.append(list(range(off, off + len(q))))
                off += len(q)
            in_group = np.zeros(n, bool)
            for g in groups:
                in_group[g] = True
            gs = np.array([len(g) for g in groups])
            lens = np.array(p["lens"])
            glens = lens[sizes_all >= MIN_GROUP]
            v = p["v"]
            Xr = np.column_stack([P[:, 0], P[:, 1], v * P[:, 2]])
            Xn = P[:, :2]
            raw = p["raw_strokes"]
            gaps = np.array([raw[i + 1][0, 2] - raw[i][-1, 2] for i in range(len(raw) - 1)])
            durs = np.array([s[-1, 2] - s[0, 2] for s in raw])
            base = f"{cat}{p['key_id']}"
            for k in KS:
                name = f"{base}_k{k}"
                E = knn_edges(Xr, k)
                N = knn_edges(Xn, k)
                for suffix, edges in (("", E), ("_null", N)):
                    od = os.path.join(args.out, name + suffix)
                    os.makedirs(od, exist_ok=True)
                    write_edges(os.path.join(od, "edges.tsv"), edges)
                    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                        for g in groups:
                            fh.write(" ".join(map(str, g)) + "\n")
                    with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                        fh.write("id\tx\ty\tt\n")
                        for i in range(n):
                            fh.write(f"{i}\t{P[i, 0]:.4f}\t{P[i, 1]:.4f}\t{P[i, 2]:.3f}\n")
                s_e = graph_stats(n, E, stroke_of, pos, in_group)
                s_n = graph_stats(n, N, stroke_of, pos, in_group)
                row = dict(graph=name, control=name + "_null", category=cat,
                           key_id=p["key_id"], countrycode=p["countrycode"],
                           timestamp=p["timestamp"], k=k, n=n, **s_e,
                           n_strokes=len(pts), n_strokes_small=int((sizes_all < MIN_GROUP).sum()),
                           n_nodes_in_small_strokes=int(sizes_all[sizes_all < MIN_GROUP].sum()),
                           n_groups=len(groups),
                           group_size_min=int(gs.min()), group_size_median=float(np.median(gs)),
                           group_size_mean=float(gs.mean()), group_size_max=int(gs.max()),
                           frac_nodes_in_group=float(in_group.mean()),
                           stroke_length_min=float(glens.min()),
                           stroke_length_median=float(np.median(glens)),
                           stroke_length_max=float(glens.max()),
                           ink_length_total=float(lens.sum()),
                           duration_ms=float(raw[-1][-1, 2]),
                           stroke_duration_ms_median=float(np.median(durs)),
                           penup_gap_ms_median=float(np.median(gaps)),
                           penup_gap_ms_min=float(gaps.min()),
                           median_speed_units_per_ms=v,
                           **{"null_" + key: val for key, val in s_n.items()})
                summary.append(row)
                print(json.dumps({x: row[x] for x in ("graph", "n", "m", "n_groups",
                      "frac_within", "null_frac_within")}))
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(args.out, "selection.json"), "w", newline="\n") as fh:
        json.dump(scan, fh, indent=1)


def cmd_table(args):
    S = json.load(open(os.path.join(args.out, "summary.json")))
    print("| graph | n | strokes (groups) | group size min-median-max | stroke length median "
          "| m | m null | mean deg | mean deg null | within | within null | consec. | consec. null |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in S:
        print(f"| {r['graph']} | {r['n']} | {r['n_strokes']} ({r['n_groups']}) | "
              f"{r['group_size_min']}-{r['group_size_median']:g}-{r['group_size_max']} | "
              f"{r['stroke_length_median']:.1f} | {r['m']} | {r['null_m']} | "
              f"{r['mean_degree']:.2f} | {r['null_mean_degree']:.2f} | "
              f"{r['frac_within']:.3f} | {r['null_frac_within']:.3f} | "
              f"{r['frac_edges_consecutive']:.3f} | {r['null_frac_edges_consecutive']:.3f} |")
    for k in KS:
        R = [r for r in S if r["k"] == k]
        for key in ("m", "mean_degree", "clustering", "frac_within", "frac_edges_consecutive",
                    "n_components"):
            a = np.array([r[key] for r in R]); b = np.array([r["null_" + key] for r in R])
            print(f"k={k} {key}: real median {np.median(a):.4g} ({a.min():.4g}-{a.max():.4g}); "
                  f"null median {np.median(b):.4g} ({b.min():.4g}-{b.max():.4g})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("build")
    a.add_argument("--raw", required=True, help="folder with <category>.part.ndjson")
    a.add_argument("--out", required=True)
    a = sub.add_parser("table")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"build": cmd_build, "table": cmd_table}[args.cmd](args)


if __name__ == "__main__":
    main()
