"""Build pixel graphs from natural images (BSDS500 test split), with the human
segmentations as ground-truth groups.

Source data (downloaded by hand with curl into its own folder, never executed):
  * BSDS500 (Arbelaez, Maire, Fowlkes, Malik 2011), official archive
      https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/BSR/BSR_bsds500.tgz
    70,763,455 bytes, sha256 97e49d31...c1a41af8e (full value in data_image/README.md).
    Only BSR/BSDS500/data/images/test/*.jpg and
    BSR/BSDS500/data/groundTruth/test/*.mat are read.

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/image_build.py build --bsr <dir that holds BSR/> --out analysis/data_image
  python -I analysis/image_build.py null  --out analysis/data_image      # only rewrite *_null/
(run from paper/; any working directory works, all paths are arguments).

Rules (fixed before any graph was built; no detection was run):
  Selection: the 200 test-split ids sorted as integers, then
    random.Random(SEED).sample(ids, N_IMAGES), then sorted again.
  Grid: centre-crop the image to a multiple of BLOCK in each direction
    (481x321 -> 480x318), cut into BLOCKxBLOCK blocks; one node per block
    (80x53 = 4240 nodes; portrait images 53x80). Node id = row * W + col.
  Colour: sRGB -> CIELAB (skimage rgb2lab, D65) per pixel, then the mean
    L, a, b over the block's 36 pixels.
  Edges (k in KS): for node i, candidates = nodes j != i with Chebyshev grid
    distance <= RADIUS. Keep the k candidates with the smallest CIELAB
    distance |lab_i - lab_j| (Delta E 1976); ties by smaller Euclidean grid
    distance, then smaller j. Undirected edge set = union over i.
  Control (<name>_null): same nodes and groups; for node i keep the k nodes
    j != i with the smallest Euclidean grid distance (no radius limit),
    ties by smaller j; union over i. Pure grid, no colour.
  Groups: for each annotator, each block gets the label held by most of its
    36 pixels (ties to the smaller label); each label gives one node set.
    All annotators' sets are pooled; exact duplicate node sets are kept once
    (first occurrence, annotator order of the .mat file), sets with < 3
    nodes are dropped, and a set holding all n nodes would be dropped.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import shutil
from collections import Counter

import numpy as np

SEED = 2026
N_IMAGES = 15
BLOCK = 6
RADIUS = 3
KS = (8, 16)
MIN_GROUP = 3


# ---------------------------------------------------------------- inputs
def test_ids(bsr):
    d = os.path.join(bsr, "BSR", "BSDS500", "data", "images", "test")
    ids = sorted((int(f[:-4]) for f in os.listdir(d) if f.endswith(".jpg")))
    gt = os.path.join(bsr, "BSR", "BSDS500", "data", "groundTruth", "test")
    gids = sorted(int(f[:-4]) for f in os.listdir(gt) if f.endswith(".mat"))
    assert ids == gids and len(ids) == 200, (len(ids), len(gids))
    return ids


def select(ids):
    return sorted(random.Random(SEED).sample(sorted(ids), N_IMAGES))


def load(bsr, iid):
    from PIL import Image
    import scipy.io as sio
    from skimage.color import rgb2lab
    p = os.path.join(bsr, "BSR", "BSDS500", "data")
    im = np.asarray(Image.open(os.path.join(p, "images", "test", f"{iid}.jpg"))
                    .convert("RGB"), dtype=np.float64) / 255.0
    mat = sio.loadmat(os.path.join(p, "groundTruth", "test", f"{iid}.mat"))
    segs = [np.asarray(mat["groundTruth"][0, a]["Segmentation"][0, 0], np.int64)
            for a in range(mat["groundTruth"].shape[1])]
    for s in segs:
        assert s.shape == im.shape[:2], (iid, s.shape, im.shape)
    return rgb2lab(im), segs


def crop_box(h, w):
    H, W = h // BLOCK, w // BLOCK
    r0, c0 = (h - H * BLOCK) // 2, (w - W * BLOCK) // 2
    return H, W, r0, c0


def block_mean(x, H, W, r0, c0):
    x = x[r0:r0 + H * BLOCK, c0:c0 + W * BLOCK]
    return x.reshape(H, BLOCK, W, BLOCK, -1).mean(axis=(1, 3))


def block_majority(seg, H, W, r0, c0):
    s = seg[r0:r0 + H * BLOCK, c0:c0 + W * BLOCK]
    s = s.reshape(H, BLOCK, W, BLOCK).transpose(0, 2, 1, 3).reshape(H, W, -1)
    out = np.empty((H, W), np.int64)
    for r in range(H):
        for c in range(W):
            vals, cnt = np.unique(s[r, c], return_counts=True)  # vals sorted
            out[r, c] = vals[np.argmax(cnt)]                     # first max = smaller label
    return out


# ---------------------------------------------------------------- edges
def grid_offsets(radius):
    return [(dr, dc) for dr in range(-radius, radius + 1)
            for dc in range(-radius, radius + 1) if (dr, dc) != (0, 0)]


def colour_edges(lab, H, W, k):
    """lab: (H*W, 3). Top-k by Delta E within Chebyshev RADIUS."""
    offs = grid_offsets(RADIUS)
    E = set()
    for r in range(H):
        for c in range(W):
            i = r * W + c
            cand = [(rr * W + cc, dr * dr + dc * dc)
                    for dr, dc in offs
                    for rr, cc in [(r + dr, c + dc)]
                    if 0 <= rr < H and 0 <= cc < W]
            j = np.array([x[0] for x in cand])
            d2 = np.array([x[1] for x in cand])
            de = np.sqrt(((lab[j] - lab[i]) ** 2).sum(1))
            order = np.lexsort((j, d2, de))          # last key is primary
            for jj in j[order[:k]]:
                jj = int(jj)
                E.add((min(i, jj), max(i, jj)))
    return sorted(E)


def null_edges(H, W, k):
    """Top-k nearest grid nodes (Euclidean), ties by smaller id; union."""
    rad = 1
    while (2 * rad + 1) ** 2 - 1 < 4 * k:            # window large enough even at a corner
        rad += 1
    rad += 1
    offs = grid_offsets(rad)
    E = set()
    for r in range(H):
        for c in range(W):
            i = r * W + c
            cand = sorted((dr * dr + dc * dc, (r + dr) * W + (c + dc))
                          for dr, dc in offs
                          if 0 <= r + dr < H and 0 <= c + dc < W)
            # a corner node needs its k nearest; check the window held them
            assert len(cand) >= k
            kth = cand[k - 1][0]
            assert kth <= rad * rad, "window too small"
            for _, j in cand[:k]:
                E.add((min(i, j), max(i, j)))
    return sorted(E)


# ---------------------------------------------------------------- groups
def build_groups(maj_list, n):
    groups, seen = [], set()
    n_small = n_dup = n_full = 0
    origin = []
    for a, m in enumerate(maj_list):
        flat = m.ravel()
        for lab in np.unique(flat):
            mem = tuple(int(x) for x in np.nonzero(flat == lab)[0])
            if len(mem) < MIN_GROUP:
                n_small += 1
                continue
            if len(mem) == n:
                n_full += 1
                continue
            if mem in seen:
                n_dup += 1
                continue
            seen.add(mem)
            groups.append(mem)
            origin.append(a)
    return groups, origin, dict(small=n_small, dup=n_dup, full=n_full)


# ---------------------------------------------------------------- stats
def graph_stats(n, edges, groups, maj_flat, H, W):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    e = np.array(edges)
    within = float(np.mean([bool(memb[a] & memb[b]) for a, b in edges]))
    per_ann = [float(np.mean(m[e[:, 0]] == m[e[:, 1]])) for m in maj_flat]
    ra, ca = np.divmod(e[:, 0], W)
    rb, cb = np.divmod(e[:, 1], W)
    cheb = np.maximum(abs(ra - rb), abs(ca - cb))
    deg = np.bincount(e.ravel(), minlength=n)
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                degree_min=int(deg.min()), degree_max=int(deg.max()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                frac_edges_cheb1=float(np.mean(cheb == 1)),
                frac_edges_within_group=within,
                frac_edges_same_segment_annot_mean=float(np.mean(per_ann)))


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


def cmd_build(args):
    ids = test_ids(args.bsr)
    chosen = select(ids)
    print("selected", chosen)
    os.makedirs(args.out, exist_ok=True)
    summary = []
    for iid in chosen:
        lab_full, segs = load(args.bsr, iid)
        h, w = lab_full.shape[:2]
        H, W, r0, c0 = crop_box(h, w)
        n = H * W
        lab = block_mean(lab_full, H, W, r0, c0).reshape(n, 3)
        maj = [block_majority(s, H, W, r0, c0) for s in segs]
        maj_flat = [m.ravel() for m in maj]
        groups, origin, dropped = build_groups(maj, n)
        sizes = [len(g) for g in groups]
        covered = set().union(*map(set, groups))
        cnt = Counter(u for g in groups for u in g)
        gsets = [set(g) for g in groups]
        near_dup = sum(1 for x in range(len(gsets)) for y in range(x + 1, len(gsets))
                       if len(gsets[x] & gsets[y]) >= 0.9 * len(gsets[x] | gsets[y]))
        seg_full = [len(np.unique(s)) for s in segs]
        seg_grid = [len(np.unique(m)) for m in maj]
        rows = np.repeat(np.arange(H), W)
        cols = np.tile(np.arange(W), H)
        for k in KS:
            name = f"bsds{iid}_k{k}"
            for suffix in ("", "_null"):
                od = os.path.join(args.out, name + suffix)
                os.makedirs(od, exist_ok=True)
                with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                    for g in groups:
                        fh.write(" ".join(map(str, g)) + "\n")
                with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                    fh.write("id\trow\tcol\tL\ta\tb\tsegment_per_annotator\n")
                    for i in range(n):
                        sl = ",".join(str(int(m[i])) for m in maj_flat)
                        fh.write(f"{i}\t{rows[i]}\t{cols[i]}\t{lab[i, 0]:.3f}\t"
                                 f"{lab[i, 1]:.3f}\t{lab[i, 2]:.3f}\t{sl}\n")
            E = colour_edges(lab, H, W, k)
            N = null_edges(H, W, k)
            write_edges(os.path.join(args.out, name, "edges.tsv"), E)
            write_edges(os.path.join(args.out, name + "_null", "edges.tsv"), N)
            s_e = graph_stats(n, E, groups, maj_flat, H, W)
            s_n = graph_stats(n, N, groups, maj_flat, H, W)
            row = dict(graph=name, control=name + "_null", image_id=iid, k=k,
                       image_h=h, image_w=w, grid_h=H, grid_w=W,
                       crop_offset_rc=[r0, c0], n=n, **s_e,
                       n_annotators=len(segs), segments_per_annotator_full=seg_full,
                       segments_per_annotator_grid=seg_grid,
                       n_groups=len(groups),
                       n_groups_dropped_small=dropped["small"],
                       n_groups_dropped_dup=dropped["dup"],
                       n_groups_dropped_full=dropped["full"],
                       groups_per_annotator=[origin.count(a) for a in range(len(segs))],
                       group_size_min=min(sizes), group_size_median=float(np.median(sizes)),
                       group_size_max=max(sizes),
                       frac_nodes_in_group=len(covered) / n,
                       groups_per_node_mean=float(np.mean([cnt[u] for u in range(n)])),
                       group_pairs_jaccard_ge_0_9=near_dup,
                       **{"null_" + key: v for key, v in s_n.items()})
            summary.append(row)
            print(json.dumps({x: row[x] for x in ("graph", "n", "m", "mean_degree",
                  "clustering", "n_groups", "frac_edges_within_group",
                  "null_frac_edges_within_group")}))
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def cmd_null(args):
    """Rewrite only the *_null folders from the real folders' nodes.tsv."""
    for name in sorted(os.listdir(args.out)):
        src = os.path.join(args.out, name)
        if name.endswith("_null") or not os.path.isfile(os.path.join(src, "edges.tsv")):
            continue
        k = int(name.rsplit("_k", 1)[1])
        rows = [l.split("\t") for l in open(os.path.join(src, "nodes.tsv")).read().splitlines()[1:]]
        r = np.array([int(x[1]) for x in rows])
        c = np.array([int(x[2]) for x in rows])
        H, W = r.max() + 1, c.max() + 1
        assert (np.array([int(x[0]) for x in rows]) == r * W + c).all()
        dst = os.path.join(args.out, name + "_null")
        os.makedirs(dst, exist_ok=True)
        N = null_edges(H, W, k)
        write_edges(os.path.join(dst, "edges.tsv"), N)
        for f in ("groups.txt", "nodes.tsv"):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
        print(name + "_null", H * W, len(N))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("build")
    a.add_argument("--bsr", required=True, help="folder that contains BSR/BSDS500/")
    a.add_argument("--out", required=True)
    a = sub.add_parser("null")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"build": cmd_build, "null": cmd_null}[args.cmd](args)


if __name__ == "__main__":
    main()
