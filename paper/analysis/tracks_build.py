"""Build detection graphs from multi-object tracking ground truth (MOT17 train
split), with the track IDs as ground-truth groups (long thin chains in time).

Source data (downloaded by hand with curl into its own folder, never executed):
  * MOT17 labels (Milan, Leal-Taixe, Reid, Roth, Schindler 2016), official archive
      https://motchallenge.net/data/MOT17Labels.zip
    10,107,022 bytes, sha256 0aa79322...dc82b2b9 (full value in data_tracks/README.md).
    Only train/MOT17-XX-FRCNN/gt/gt.txt and train/MOT17-XX-FRCNN/seqinfo.ini
    are read (the DPM / FRCNN / SDP copies of each gt.txt are byte-identical).
    No image is used.

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/tracks_build.py build --mot <dir that holds train/> --out analysis/data_tracks
(run from paper/; any working directory works, all paths are arguments).

Rules (fixed before any graph was built; no detection was run):
  Nodes: rows of gt.txt with flag (col 7) == 1, class (col 8) == 1
    (pedestrian) and visibility (col 9) >= MIN_VIS. One node per box.
  Windows: per sequence, N = kept boxes, c = floor(N / TARGET + 0.5)
    (sequence skipped if c == 0). The frames 1..seqLength are cut into c runs
    of consecutive frames: cut j (j = 1..c-1) goes after the first frame at
    which the running box count reaches >= j * N / c.
  Selection: windows listed as (sequence, j), sequences sorted. With
    rng = random.Random(SEED): first one window per sequence,
    rng.choice(windows of that sequence), sequences in sorted order; then
    rng.sample(remaining windows in sorted order, N_WINDOWS - picked).
  Node order inside a window: sorted by (frame, x, y, w, h, line number in
    gt.txt); id = rank, 0-based. The track id plays no part in the order.
  Candidates of node a: every other node b (b != a) with |frame_b - frame_a|
    <= DT (same frame included) and centre distance |c_a - c_b| <= RAD * h_a
    (h_a = height of a's box; centre = (x + w/2, y + h/2)).
  Edges (k in KS): node a keeps the k candidates with the largest IoU of the
    two boxes as given (no shift), ties by smaller centre distance, then by
    smaller id (all candidates if fewer than k). Undirected edge set = union.
  Control (<name>_null): same candidates; keep the k with the smallest centre
    distance, ties by smaller id; union. Same nodes.tsv and groups.txt.
  Groups: track ids with >= MIN_GROUP nodes in the window, one group each
    (disjoint). Nodes of shorter tracks stay in the graph but in no group.
"""
from __future__ import annotations

import argparse
import configparser
import json
import os
import random

import numpy as np

SEED = 2026
TARGET = 3500
N_WINDOWS = 15
MIN_VIS = 0.25
DT = 3
RAD = 1.5
KS = (8, 16)
MIN_GROUP = 3


# ---------------------------------------------------------------- inputs
def sequences(mot):
    d = os.path.join(mot, "train")
    seqs = sorted(s for s in os.listdir(d) if s.endswith("-FRCNN"))
    assert len(seqs) == 7, seqs
    return seqs


def load(mot, seq):
    d = os.path.join(mot, "train", seq)
    cp = configparser.ConfigParser()
    cp.read(os.path.join(d, "seqinfo.ini"))
    info = dict(seq_length=int(cp["Sequence"]["seqLength"]),
                im_w=int(cp["Sequence"]["imWidth"]),
                im_h=int(cp["Sequence"]["imHeight"]),
                frame_rate=int(cp["Sequence"]["frameRate"]))
    raw = np.loadtxt(os.path.join(d, "gt", "gt.txt"), delimiter=",", ndmin=2)
    assert raw.shape[1] == 9, raw.shape
    keep = (raw[:, 6] == 1) & (raw[:, 7] == 1) & (raw[:, 8] >= MIN_VIS)
    line = np.nonzero(keep)[0]
    r = raw[keep]
    det = dict(frame=r[:, 0].astype(np.int64), tid=r[:, 1].astype(np.int64),
               x=r[:, 2], y=r[:, 3], w=r[:, 4], h=r[:, 5], line=line)
    assert det["frame"].min() >= 1 and det["frame"].max() <= info["seq_length"]
    info.update(n_rows=len(raw), n_kept=int(keep.sum()))
    return det, info


def windows(frame, seq_length):
    N = len(frame)
    c = int(np.floor(N / TARGET + 0.5))
    if c == 0:
        return []
    cnt = np.bincount(frame, minlength=seq_length + 1)[1:]
    cum = np.cumsum(cnt)
    cuts = [0]
    for j in range(1, c):
        f = int(np.nonzero(cum >= j * N / c)[0][0]) + 1      # 1-based frame
        cuts.append(f)
    cuts.append(seq_length)
    return [(cuts[j] + 1, cuts[j + 1]) for j in range(c)]   # inclusive frames


def select(win):
    """win: dict seq -> list of windows. Returns sorted list of (seq, j)."""
    rng = random.Random(SEED)
    picked = []
    for s in sorted(win):
        if win[s]:
            picked.append((s, rng.choice(range(len(win[s])))))
    rest = sorted((s, j) for s in win for j in range(len(win[s])) if (s, j) not in picked)
    picked += rng.sample(rest, min(N_WINDOWS - len(picked), len(rest)))
    return sorted(picked)


# ---------------------------------------------------------------- edges
def build_edges(fr, x, y, w, h, ks, gid):
    """Returns {k: (real_edges, null_edges)}, candidates per node, and the
    unordered same-group candidate pairs. Nodes sorted by frame."""
    n = len(fr)
    cx, cy = x + w / 2, y + h / 2
    x2, y2 = x + w, y + h
    area = w * h
    ids = np.arange(n)
    real = {k: set() for k in ks}
    null = {k: set() for k in ks}
    ncand = np.zeros(n, np.int64)
    gpairs = set()
    for f in np.unique(fr):
        A = np.nonzero(fr == f)[0]
        B = np.nonzero(np.abs(fr - f) <= DT)[0]
        d = np.hypot(cx[A, None] - cx[None, B], cy[A, None] - cy[None, B])
        iw = np.clip(np.minimum(x2[A, None], x2[None, B]) - np.maximum(x[A, None], x[None, B]), 0, None)
        ih = np.clip(np.minimum(y2[A, None], y2[None, B]) - np.maximum(y[A, None], y[None, B]), 0, None)
        inter = iw * ih
        iou = inter / (area[A, None] + area[None, B] - inter)
        for ia, a in enumerate(A):
            m = (d[ia] <= RAD * h[a]) & (B != a)
            j, dj, uj = B[m], d[ia][m], iou[ia][m]
            ncand[a] = len(j)
            if gid[a] >= 0:
                for b in j[gid[j] == gid[a]]:
                    gpairs.add((min(a, int(b)), max(a, int(b))))
            o_real = np.lexsort((j, dj, -uj))                 # last key primary
            o_null = np.lexsort((j, dj))
            for k in ks:
                for b in j[o_real[:k]]:
                    real[k].add((min(a, int(b)), max(a, int(b))))
                for b in j[o_null[:k]]:
                    null[k].add((min(a, int(b)), max(a, int(b))))
    return {k: (sorted(real[k]), sorted(null[k])) for k in ks}, ncand, gpairs


# ---------------------------------------------------------------- stats
def graph_stats(n, edges, fr, gid, groups, gpairs):
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    e = np.array(edges, dtype=np.int64).reshape(-1, 2)
    A = sp.coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    A = ((A + A.T) > 0).astype(np.float64).tocsr()
    deg = np.asarray(A.sum(1)).ravel()
    tri = np.asarray((A @ A).multiply(A).sum(1)).ravel() / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        cc = np.where(deg >= 2, 2 * tri / (deg * (deg - 1)), 0.0)
    ncomp = connected_components(A, directed=False)[0]
    same = (gid[e[:, 0]] >= 0) & (gid[e[:, 0]] == gid[e[:, 1]])
    dt = np.abs(fr[e[:, 0]] - fr[e[:, 1]])
    pieces = np.array([connected_components(A[g][:, g], directed=False)[0]
                       for g in groups])
    linked = len(gpairs & set(map(tuple, e.tolist())))
    return dict(m=len(e), mean_degree=2 * len(e) / n,
                degree_min=int(deg.min()), degree_max=int(deg.max()),
                clustering=float(cc.mean()), n_components=int(ncomp),
                frac_edges_within_group=float(same.mean()),
                frac_edges_same_frame=float(np.mean(dt == 0)),
                mean_edge_dt=float(dt.mean()),
                frac_groups_connected=float(np.mean(pieces == 1)),
                group_pieces_mean=float(pieces.mean()),
                within_ceiling=float(min(1.0, len(gpairs) / len(e))),
                frac_within_candidate_pairs_linked=linked / len(gpairs))


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


def fmt(v):
    return str(int(v)) if float(v).is_integer() else f"{v:.2f}"


def cmd_build(args):
    seqs = sequences(args.mot)
    data, win = {}, {}
    for s in seqs:
        det, info = load(args.mot, s)
        data[s] = (det, info)
        win[s] = windows(det["frame"], info["seq_length"])
        print(s, info, "windows", win[s])
    chosen = select(win)
    print("selected", chosen)
    os.makedirs(args.out, exist_ok=True)
    summary = []
    for s, j in chosen:
        det, info = data[s]
        f0, f1 = win[s][j]
        sel = np.nonzero((det["frame"] >= f0) & (det["frame"] <= f1))[0]
        order = np.lexsort((det["line"][sel], det["h"][sel], det["w"][sel],
                            det["y"][sel], det["x"][sel], det["frame"][sel]))
        sel = sel[order]
        fr, tid = det["frame"][sel], det["tid"][sel]
        x, y, w, h = det["x"][sel], det["y"][sel], det["w"][sel], det["h"][sel]
        n = len(sel)
        # one box per track per frame
        assert len(set(zip(fr.tolist(), tid.tolist()))) == n, s
        groups, gid = [], -np.ones(n, np.int64)
        tracks = np.unique(tid)
        n_small = 0
        span = []
        for t in tracks:
            mem = np.nonzero(tid == t)[0]
            if len(mem) < MIN_GROUP:
                n_small += 1
                continue
            gid[mem] = len(groups)
            groups.append(mem)
            span.append(int(fr[mem].max() - fr[mem].min() + 1))
        sizes = np.array([len(g) for g in groups])
        res, ncand, gpairs = build_edges(fr, x, y, w, h, KS, gid)
        seqnum = s.split("-")[1]
        for k in KS:
            name = f"mot17_{seqnum}_w{j + 1}_k{k}"
            for suffix in ("", "_null"):
                od = os.path.join(args.out, name + suffix)
                os.makedirs(od, exist_ok=True)
                with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                    for g in groups:
                        fh.write(" ".join(map(str, g.tolist())) + "\n")
                with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                    fh.write("id\tframe\tx\ty\tw\th\n")
                    for i in range(n):
                        fh.write(f"{i}\t{fr[i]}\t{fmt(x[i])}\t{fmt(y[i])}\t"
                                 f"{fmt(w[i])}\t{fmt(h[i])}\n")
            E, N = res[k]
            write_edges(os.path.join(args.out, name, "edges.tsv"), E)
            write_edges(os.path.join(args.out, name + "_null", "edges.tsv"), N)
            s_e = graph_stats(n, E, fr, gid, groups, gpairs)
            s_n = graph_stats(n, N, fr, gid, groups, gpairs)
            row = dict(graph=name, control=name + "_null", sequence=s.rsplit("-", 1)[0],
                       window=j + 1, n_windows_in_sequence=len(win[s]),
                       frame_first=int(f0), frame_last=int(f1), n_frames=int(f1 - f0 + 1),
                       frame_rate=info["frame_rate"], image_w=info["im_w"],
                       image_h=info["im_h"], k=k, n=n, **s_e,
                       candidates_mean=float(ncand.mean()),
                       candidates_min=int(ncand.min()),
                       frac_nodes_lt_k_candidates=float(np.mean(ncand < k)),
                       n_within_candidate_pairs=len(gpairs),
                       n_tracks=len(tracks), n_groups=len(groups),
                       n_tracks_dropped_small=n_small,
                       group_size_min=int(sizes.min()),
                       group_size_mean=float(sizes.mean()),
                       group_size_median=float(np.median(sizes)),
                       group_size_max=int(sizes.max()),
                       group_span_frames_mean=float(np.mean(span)),
                       group_span_frames_max=int(max(span)),
                       frac_nodes_in_group=float(np.mean(gid >= 0)),
                       **{"null_" + key: v for key, v in s_n.items()})
            summary.append(row)
            print(json.dumps({x_: row[x_] for x_ in ("graph", "n", "m", "mean_degree",
                  "clustering", "n_groups", "group_size_mean",
                  "frac_edges_within_group", "null_frac_edges_within_group")}))
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("build")
    a.add_argument("--mot", required=True, help="folder that contains train/MOT17-XX-FRCNN/")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"build": cmd_build}[args.cmd](args)


if __name__ == "__main__":
    main()
