"""Build hit graphs from TrackML simulated collision events, with the particles
(tracks) as ground-truth groups: long thin chains of hits.

Source data (fetched with HTTP Range requests into its own folder, never
executed; only hits/truth CSVs of the selected events were read):
  * TrackML Particle Tracking Challenge, Zenodo record 4730167
    (doi:10.5281/zenodo.4730167, CC BY 4.0), file train_1.zip
      https://zenodo.org/records/4730167
    Only train_1/event<ID>-hits.csv and train_1/event<ID>-truth.csv of the
    selected events are read (CRC-checked by zipfile on extraction).

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/trackml_build.py ids    --out analysis/data_trackml/train_1_event_ids.txt
  python -I analysis/trackml_build.py select --ids analysis/data_trackml/train_1_event_ids.txt
  python -I analysis/trackml_build.py fetch  --ids analysis/data_trackml/train_1_event_ids.txt \
      --dest <new empty folder>
  python -I analysis/trackml_build.py build  --raw <dir with event*-hits.csv, event*-truth.csv> \
      --ids analysis/data_trackml/train_1_event_ids.txt --out analysis/data_trackml
`ids` reads only the central directory of the remote train_1.zip; `fetch`
pulls only the 16 selected members (about 71 MB) with HTTP Range requests,
standard library only, CRC-checked; nothing downloaded is executed.
(run from paper/; any working directory works, all paths are arguments).

Rules (fixed before any graph was built; no detection was run):
  Selection: the 1770 event ids of train_1.zip (train_1_event_ids.txt, read
    from the zip's central directory), sorted; rng = random.Random(SEED);
    events = sorted(rng.sample(ids, N_EVENTS)); then for each chosen event in
    that order: sector = rng.randrange(N_PHI), eta_bin = rng.randrange(2).
  Slice: barrel volumes 8, 13, 17 only (all their layers). For a hit with
    r = hypot(x, y): phi = atan2(y, x), eta = asinh(z / r) (position, not
    truth). Keep hits with -pi + s*2pi/N_PHI <= phi < -pi + (s+1)*2pi/N_PHI
    and eta in [-1, 0) (eta_bin 0) or [0, 1) (eta_bin 1). Noise hits
    (particle_id 0) stay as nodes. Node id = rank of hit_id in the slice.
  Layer order: the 10 barrel layers sorted by radius get a global index
    L = 0..9: (8,2) (8,4) (8,6) (8,8) (13,2) (13,4) (13,6) (13,8) (17,2) (17,4).
    "Adjacent" is taken in this barrel-wide order, so a track can be linked
    across the pixel / short-strip / long-strip volume boundaries (within one
    volume only, every track of >= 5 layers would be cut into pieces).
  Candidates of hit a: hits b in the slice with |L_b - L_a| in {1, 2},
    |dphi| <= DPHI (wrapped) and |z_b - z_a| <= DZ.
  Edges (k in KS): keep the k candidates with the smallest angle between
    the line through a and b and the line through the beam spot (0,0,0) and
    a (unsigned, 0..90 deg); ties by smaller 3D distance, then smaller id.
    Fewer than k candidates -> keep all. Undirected edge set = union over a.
  Control (<name>_null): same nodes, groups and candidate sets; keep the k
    candidates with the smallest 3D distance |b - a|, ties by smaller id;
    union over a.
  Groups: particles (particle_id != 0) with >= MIN_GROUP hits in the slice;
    one group per particle (a partition of the non-noise hits it covers).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import shutil

import numpy as np

ZIP_URL = "https://zenodo.org/api/records/4730167/files/train_1.zip/content"
SEED = 2026
N_EVENTS = 8
N_PHI = 8
ETA_BINS = ((-1.0, 0.0), (0.0, 1.0))
DPHI = 0.2          # rad
DZ = 500.0          # mm
LAYER_STEPS = (1, 2)
KS = (8, 16)
MIN_GROUP = 3
BARREL = [(8, 2), (8, 4), (8, 6), (8, 8), (13, 2), (13, 4), (13, 6), (13, 8), (17, 2), (17, 4)]
GLOBAL_LAYER = {vl: i for i, vl in enumerate(BARREL)}


# ---------------------------------------------------------------- remote zip
class HTTPRangeFile:
    """Read-only seekable view of a remote file, by HTTP Range requests."""

    def __init__(self, url, chunk=1 << 20):
        import urllib.request
        self._ur = urllib.request
        self.url, self.chunk = url, chunk
        with self._open(0, 0) as r:
            self.size = int(r.headers["Content-Range"].split("/")[1])
        self.pos, self._b0, self._buf = 0, 0, b""

    def _open(self, a, b):
        req = self._ur.Request(self.url, headers={"User-Agent": "python-urllib",
                                                  "Range": f"bytes={a}-{b}"})
        r = self._ur.urlopen(req, timeout=60)
        assert r.status == 206, r.status
        return r

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        if not (self._b0 <= self.pos and self.pos + n <= self._b0 + len(self._buf)):
            end = min(self.pos + max(n, self.chunk), self.size) - 1
            with self._open(self.pos, end) as r:
                self._buf = r.read()
            self._b0 = self.pos
        s = self.pos - self._b0
        out = self._buf[s:s + n]
        self.pos += len(out)
        return out


def remote_zip():
    import zipfile
    return zipfile.ZipFile(HTTPRangeFile(ZIP_URL))


def member(ev, kind):
    return f"train_1/event{ev:09d}-{kind}.csv"


# ---------------------------------------------------------------- selection
def read_ids(path):
    ids = sorted(int(x) for x in open(path).read().split())
    assert len(ids) == len(set(ids)) == 1770, len(ids)
    return ids


def select(ids):
    rng = random.Random(SEED)
    events = sorted(rng.sample(ids, N_EVENTS))
    return [(e, rng.randrange(N_PHI), rng.randrange(len(ETA_BINS))) for e in events]


def graph_name(ev, sec, eb, k):
    return f"trackml{ev}_p{sec}e{eb}_k{k}"


# ---------------------------------------------------------------- inputs
def load_slice(raw, ev, sec, eb):
    import pandas as pd
    hits = pd.read_csv(os.path.join(raw, f"event{ev:09d}-hits.csv"))
    truth = pd.read_csv(os.path.join(raw, f"event{ev:09d}-truth.csv"),
                        usecols=["hit_id", "particle_id"])
    assert len(hits) == len(truth)
    df = hits.merge(truth, on="hit_id", how="left", validate="1:1")
    assert df["particle_id"].notna().all()
    n_event = len(df)
    key = list(zip(df["volume_id"], df["layer_id"]))
    df["L"] = [GLOBAL_LAYER.get(t, -1) for t in key]
    barrel = df["volume_id"].isin([8, 13, 17]).to_numpy()
    assert (df["L"].to_numpy()[barrel] >= 0).all(), "unexpected barrel layer id"
    r = np.hypot(df["x"].to_numpy(), df["y"].to_numpy())
    phi = np.arctan2(df["y"].to_numpy(), df["x"].to_numpy())
    eta = np.arcsinh(df["z"].to_numpy() / r)
    w = 2 * math.pi / N_PHI
    lo, hi = -math.pi + sec * w, -math.pi + (sec + 1) * w
    e0, e1 = ETA_BINS[eb]
    keep = barrel & (phi >= lo) & (phi < hi) & (eta >= e0) & (eta < e1)
    s = df[keep].sort_values("hit_id").reset_index(drop=True)
    return s, n_event, int(barrel.sum())


# ---------------------------------------------------------------- edges
def candidates(xyz, L):
    """For every hit a: array of candidate ids b (sorted), per the window rule."""
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    z = xyz[:, 2]
    by_layer = {l: np.nonzero(L == l)[0] for l in range(len(BARREL))}
    out = []
    for a in range(len(L)):
        cs = []
        for st in LAYER_STEPS:
            for lb in (L[a] - st, L[a] + st):
                if lb not in by_layer:
                    continue
                b = by_layer[lb]
                dphi = np.abs((phi[b] - phi[a] + math.pi) % (2 * math.pi) - math.pi)
                cs.append(b[(dphi <= DPHI) & (np.abs(z[b] - z[a]) <= DZ)])
        out.append(np.sort(np.concatenate(cs)) if cs else np.zeros(0, np.int64))
    return out


def topk_edges(xyz, cand, k, score):
    E = set()
    for a, b in enumerate(cand):
        if len(b) == 0:
            continue
        u = xyz[b] - xyz[a]
        dist = np.sqrt((u ** 2).sum(1))
        if score == "origin":
            v = xyz[a]
            cosang = np.abs(u @ v) / (dist * np.sqrt(v @ v))
            ang = np.arccos(np.clip(cosang, 0.0, 1.0))
            order = np.lexsort((b, dist, ang))       # last key is primary
        else:
            order = np.lexsort((b, dist))
        for bb in b[order[:k]]:
            bb = int(bb)
            E.add((min(a, bb), max(a, bb)))
    return sorted(E)


# ---------------------------------------------------------------- stats
def graph_stats(n, edges, groups, pid):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    e = np.array(edges)
    within = float(np.mean((pid[e[:, 0]] == pid[e[:, 1]]) & (pid[e[:, 0]] != 0)))
    noise_end = float(np.mean((pid[e[:, 0]] == 0) | (pid[e[:, 1]] == 0)))
    deg = np.bincount(e.ravel(), minlength=n)
    conn = [nx.is_connected(G.subgraph(g)) for g in groups]
    # recall of true track segments: pairs of hits of one group on
    # consecutive global layers (|dL| == 1)
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                degree_min=int(deg.min()), degree_max=int(deg.max()),
                n_isolated=int((deg == 0).sum()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                frac_edges_within_group=within,
                frac_edges_with_noise_end=noise_end,
                frac_groups_connected=float(np.mean(conn)))


def seg_recall(edges, groups, L):
    E = set(edges)
    tot = hit = 0
    for g in groups:
        for i, a in enumerate(g):
            for b in g[i + 1:]:
                if abs(int(L[a]) - int(L[b])) == 1:
                    tot += 1
                    hit += (a, b) in E
    return hit / tot if tot else float("nan")


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


# ---------------------------------------------------------------- commands
def cmd_ids(args):
    import re
    names = remote_zip().namelist()
    ids = sorted({int(m.group(1)) for x in names
                  for m in [re.fullmatch(r"train_1/event(\d{9})-hits\.csv", x)] if m})
    with open(args.out, "w", newline="\n") as fh:
        fh.write("\n".join(map(str, ids)) + "\n")
    print(len(ids), "event ids ->", args.out)


def cmd_select(args):
    for ev, sec, eb in select(read_ids(args.ids)):
        print(ev, sec, eb, member(ev, "hits"), member(ev, "truth"))


def cmd_fetch(args):
    os.makedirs(args.dest, exist_ok=True)
    z = remote_zip()
    for ev, _, _ in select(read_ids(args.ids)):
        for kind in ("hits", "truth"):
            info = z.getinfo(member(ev, kind))
            dst = os.path.join(args.dest, os.path.basename(info.filename))
            assert not os.path.exists(dst), dst
            with z.open(info) as src, open(dst, "wb") as fh:   # zipfile checks the CRC
                shutil.copyfileobj(src, fh, 1 << 20)
            print(info.filename, info.file_size, f"crc32={info.CRC:08x}", flush=True)


def cmd_build(args):
    chosen = select(read_ids(args.ids))
    print("selected", chosen)
    os.makedirs(args.out, exist_ok=True)
    summary = []
    for ev, sec, eb in chosen:
        s, n_event, n_barrel = load_slice(args.raw, ev, sec, eb)
        n = len(s)
        xyz = s[["x", "y", "z"]].to_numpy(np.float64)
        L = s["L"].to_numpy()
        pid = s["particle_id"].to_numpy(np.int64)
        cnt = {}
        for i, p in enumerate(pid):
            if p != 0:
                cnt.setdefault(int(p), []).append(i)
        groups = sorted((tuple(v) for v in cnt.values() if len(v) >= MIN_GROUP),
                        key=lambda g: g[0])
        sizes = [len(g) for g in groups]
        covered = sum(sizes)
        n_noise = int((pid == 0).sum())
        n_small_hits = sum(len(v) for v in cnt.values() if len(v) < MIN_GROUP)
        layers_per_group = [len(set(L[list(g)].tolist())) for g in groups]
        cand = candidates(xyz, L)
        ncand = np.array([len(c) for c in cand])
        for k in KS:
            name = graph_name(ev, sec, eb, k)
            for suffix in ("", "_null"):
                od = os.path.join(args.out, name + suffix)
                os.makedirs(od, exist_ok=True)
                with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                    for g in groups:
                        fh.write(" ".join(map(str, g)) + "\n")
                with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                    fh.write("id\tx\ty\tz\tvolume\tlayer\thit_id\n")
                    for i in range(n):
                        fh.write(f"{i}\t{xyz[i, 0]:.4f}\t{xyz[i, 1]:.4f}\t{xyz[i, 2]:.4f}\t"
                                 f"{int(s['volume_id'][i])}\t{int(s['layer_id'][i])}\t"
                                 f"{int(s['hit_id'][i])}\n")
            E = topk_edges(xyz, cand, k, "origin")
            N = topk_edges(xyz, cand, k, "dist")
            write_edges(os.path.join(args.out, name, "edges.tsv"), E)
            write_edges(os.path.join(args.out, name + "_null", "edges.tsv"), N)
            s_e = graph_stats(n, E, groups, pid)
            s_n = graph_stats(n, N, groups, pid)
            s_e["frac_track_segments_recovered"] = seg_recall(E, groups, L)
            s_n["frac_track_segments_recovered"] = seg_recall(N, groups, L)
            row = dict(graph=name, control=name + "_null", event=ev, phi_sector=sec,
                       eta_bin=eb, eta_range=list(ETA_BINS[eb]),
                       phi_range=[-math.pi + sec * 2 * math.pi / N_PHI,
                                  -math.pi + (sec + 1) * 2 * math.pi / N_PHI],
                       k=k, n_hits_event=n_event, n_hits_barrel=n_barrel, n=n, **s_e,
                       n_noise_hits=n_noise,
                       n_hits_in_particles_lt3=n_small_hits,
                       n_particles_in_slice=len(cnt),
                       n_groups=len(groups),
                       group_size_min=min(sizes), group_size_median=float(np.median(sizes)),
                       group_size_mean=float(np.mean(sizes)), group_size_max=max(sizes),
                       layers_per_group_mean=float(np.mean(layers_per_group)),
                       frac_nodes_in_group=covered / n,
                       candidates_per_node_min=int(ncand.min()),
                       candidates_per_node_median=float(np.median(ncand)),
                       n_nodes_fewer_than_k_candidates=int((ncand < k).sum()),
                       **{"null_" + key: v for key, v in s_n.items()})
            summary.append(row)
            print(json.dumps({x: row[x] for x in ("graph", "n", "m", "mean_degree",
                  "n_groups", "group_size_mean", "frac_edges_within_group",
                  "null_frac_edges_within_group")}), flush=True)
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("ids")
    a.add_argument("--out", required=True)
    a = sub.add_parser("select")
    a.add_argument("--ids", required=True)
    a = sub.add_parser("fetch")
    a.add_argument("--ids", required=True)
    a.add_argument("--dest", required=True, help="new empty folder for the CSVs")
    a = sub.add_parser("build")
    a.add_argument("--raw", required=True, help="folder with event*-hits.csv and -truth.csv")
    a.add_argument("--ids", required=True)
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"ids": cmd_ids, "select": cmd_select, "fetch": cmd_fetch,
     "build": cmd_build}[args.cmd](args)


if __name__ == "__main__":
    main()
