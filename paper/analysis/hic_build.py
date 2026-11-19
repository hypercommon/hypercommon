"""Build small real networks from Hi-C contact maps, with Arrowhead TADs as groups.

Source data (read remotely, never executed):
  * Rao et al. 2014 (Cell 159:1665), GM12878 in situ Hi-C, primary+replicate,
    MAPQ>=30 Juicer .hic file on GEO GSE63525:
      HIC_URL below (40 GB; only the byte ranges we need are fetched with
      HTTP Range requests).
  * Arrowhead contact-domain list from the same GEO series (ARROWHEAD_URL).

The .hic (Juicer format v7) reader is written here from the published format
spec (github.com/aidenlab/hic-format) because the hicstraw C++ package does
not build on this Windows machine (needs libcurl headers).

Usage (run with an isolated interpreter, e.g. `python -I hic_build.py ...`):
  python -I hic_build.py probe
  python -I hic_build.py fetch  --cache D:/tmp/hic_cache --chroms 17 18 19 20 21 22
  python -I hic_build.py build  --cache D:/tmp/hic_cache --domains <arrowhead.txt.gz>
                                --out <paper/analysis/data_hic> --chroms 17 ... 22

Edge rule (per chromosome, 25 kb bins, KR-normalised observed contacts):
  1. obs[i,j] = raw[i,j] / (kr[i] * kr[j]); bins with NaN/zero KR factor or an
     all-zero raw row are dropped ("no coverage"); the rest are renumbered
     0..n-1 in genomic order.
  2. For each retained bin i, keep the k retained bins j != i with the largest
     obs[i,j] (only obs > 0; ties broken by smaller genomic distance, then by
     smaller j). Undirected edge set = union over i. Unweighted.
  Diagnostics only (not used for edges): the same top-k rule on
  observed/expected (obs / mean obs at that distance), and on a pure
  genomic-distance "band" null (score = 1/(1+|i-j|)), are summarised per graph.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import struct
import sys
import time
import urllib.request
import zlib

import numpy as np

HIC_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/"
           "GSE63525_GM12878_insitu_primary%2Breplicate_combined_30.hic")
ARROWHEAD_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/"
                 "GSE63525_GM12878_primary%2Breplicate_Arrowhead_domainlist.txt.gz")
UA = "hic-graph-builder/0.1 (python urllib)"


# ----------------------------------------------------------------------------
# Remote random-access file
# ----------------------------------------------------------------------------
class RemoteFile:
    """Seekable read-only view of an HTTP resource via Range requests."""

    CHUNK = 1 << 18  # 256 KiB cache granularity for small reads

    def __init__(self, url):
        self.url = url
        self.pos = 0
        self.cache = {}
        self.bytes_fetched = 0
        self.n_requests = 0

    def fetch(self, start, end):  # inclusive start, exclusive end
        for attempt in range(6):
            try:
                req = urllib.request.Request(
                    self.url, headers={"Range": f"bytes={start}-{end - 1}",
                                       "User-Agent": UA})
                with urllib.request.urlopen(req, timeout=120) as r:
                    if r.status != 206:
                        raise IOError(f"expected 206, got {r.status}")
                    data = r.read()
                if len(data) != end - start:
                    raise IOError(f"short read {len(data)} != {end - start}")
                self.bytes_fetched += len(data)
                self.n_requests += 1
                return data
            except Exception as e:  # noqa: BLE001
                wait = 2 ** attempt
                print(f"  fetch {start}-{end} failed ({e}); retry in {wait}s",
                      file=sys.stderr)
                time.sleep(wait)
        raise IOError(f"giving up on range {start}-{end}")

    def _chunk(self, idx):
        if idx not in self.cache:
            if len(self.cache) > 256:
                self.cache.clear()
            s = idx * self.CHUNK
            self.cache[idx] = self.fetch(s, s + self.CHUNK)
        return self.cache[idx]

    def seek(self, pos, whence=0):
        self.pos = pos if whence == 0 else self.pos + pos

    def read(self, n):
        out = bytearray()
        while n > 0:
            idx, off = divmod(self.pos, self.CHUNK)
            c = self._chunk(idx)
            take = c[off:off + n]
            if not take:
                raise EOFError
            out += take
            self.pos += len(take)
            n -= len(take)
        return bytes(out)

    # typed helpers
    def i32(self):
        return struct.unpack("<i", self.read(4))[0]

    def i64(self):
        return struct.unpack("<q", self.read(8))[0]

    def f32(self):
        return struct.unpack("<f", self.read(4))[0]

    def cstr(self):
        buf = bytearray()
        while True:
            b = self.read(1)
            if b == b"\0":
                return buf.decode("utf-8")
            buf += b


# ----------------------------------------------------------------------------
# .hic v7/v8 reader (only what we need)
# ----------------------------------------------------------------------------
class HicFile:
    def __init__(self, url):
        self.f = RemoteFile(url)
        self._read_header()
        self._read_footer()

    def _read_header(self):
        f = self.f
        f.seek(0)
        if f.read(4) != b"HIC\0":
            raise ValueError("not a .hic file")
        self.version = f.i32()
        if self.version not in (7, 8):
            raise ValueError(f"unsupported .hic version {self.version}")
        self.master = f.i64()
        self.genome = f.cstr()
        n_attr = f.i32()
        self.attributes = {}
        for _ in range(n_attr):
            k = f.cstr()
            v = f.cstr()
            self.attributes[k] = v
        n_chr = f.i32()
        self.chroms = []
        for _ in range(n_chr):
            name = f.cstr()
            length = f.i32()
            self.chroms.append((name, length))
        n_bp = f.i32()
        self.bp_res = [f.i32() for _ in range(n_bp)]
        n_fr = f.i32()
        self.frag_res = [f.i32() for _ in range(n_fr)]

    def _skip_doubles(self, n):
        self.f.seek(8 * n, 1)

    def _read_footer(self):
        f = self.f
        f.seek(self.master)
        f.i32()  # nBytes
        n = f.i32()
        self.matrix_index = {}
        for _ in range(n):
            key = f.cstr()
            pos = f.i64()
            size = f.i32()
            self.matrix_index[key] = (pos, size)
        # expected values (no normalisation)
        self.expected = {}
        n_exp = f.i32()
        for _ in range(n_exp):
            unit = f.cstr()
            bs = f.i32()
            nv = f.i32()
            if unit == "BP" and bs in (25000, 50000):
                self.expected[("NONE", unit, bs)] = np.frombuffer(
                    f.read(8 * nv), "<f8").copy()
            else:
                self._skip_doubles(nv)
            nf = f.i32()
            f.seek(12 * nf, 1)
        # normalised expected values
        n_exp = f.i32()
        for _ in range(n_exp):
            typ = f.cstr()
            unit = f.cstr()
            bs = f.i32()
            nv = f.i32()
            if unit == "BP" and bs in (25000, 50000):
                self.expected[(typ, unit, bs)] = np.frombuffer(
                    f.read(8 * nv), "<f8").copy()
            else:
                self._skip_doubles(nv)
            nf = f.i32()
            f.seek(12 * nf, 1)
        # normalisation vector index
        n = f.i32()
        self.norm_index = {}
        for _ in range(n):
            typ = f.cstr()
            ci = f.i32()
            unit = f.cstr()
            bs = f.i32()
            pos = f.i64()
            size = f.i32()
            self.norm_index[(typ, ci, unit, bs)] = (pos, size)

    def chrom_index(self, name):
        for i, (n, _) in enumerate(self.chroms):
            if n == name:
                return i
        raise KeyError(name)

    def norm_vector(self, norm, ci, bs):
        pos, size = self.norm_index[(norm, ci, "BP", bs)]
        data = self.f.fetch(pos, pos + size)
        nv = struct.unpack("<i", data[:4])[0]
        return np.frombuffer(data[4:4 + 8 * nv], "<f8").copy()

    def intra_matrix(self, ci, bs):
        """Return raw counts (x, y, v) for all upper-triangle records."""
        f = self.f
        key = f"{ci}_{ci}"
        pos, _ = self.matrix_index[key]
        f.seek(pos)
        f.i32()
        f.i32()
        n_res = f.i32()
        blocks = None
        for _ in range(n_res):
            unit = f.cstr()
            f.i32()  # zoom index
            f.read(16)  # sumCounts, occupiedCellCount, stdDev, percent95
            binsize = f.i32()
            f.i32()  # blockBinCount
            f.i32()  # blockColumnCount
            nb = f.i32()
            raw = f.read(16 * nb)
            if unit == "BP" and binsize == bs:
                arr = np.frombuffer(raw, dtype=[("num", "<i4"), ("pos", "<i8"),
                                                ("size", "<i4")])
                blocks = arr.copy()
                break
        if blocks is None:
            raise KeyError(f"no BP {bs} zoom for chrom {ci}")
        order = np.argsort(blocks["pos"])
        blocks = blocks[order]
        xs, ys, vs = [], [], []
        # fetch contiguous runs of blocks in <= 32 MB requests
        i = 0
        nb = len(blocks)
        while i < nb:
            start = int(blocks["pos"][i])
            j = i
            end = start + int(blocks["size"][i])
            while (j + 1 < nb and int(blocks["pos"][j + 1]) - end < 65536
                   and int(blocks["pos"][j + 1]) + int(blocks["size"][j + 1])
                   - start <= 32 << 20):
                j += 1
                end = max(end, int(blocks["pos"][j]) + int(blocks["size"][j]))
            buf = f.fetch(start, end)
            for b in range(i, j + 1):
                o = int(blocks["pos"][b]) - start
                comp = buf[o:o + int(blocks["size"][b])]
                x, y, v = decode_block(zlib.decompress(comp), self.version)
                xs.append(x)
                ys.append(y)
                vs.append(v)
            i = j + 1
        return (np.concatenate(xs), np.concatenate(ys), np.concatenate(vs))


def decode_block(u, version):
    n_rec = struct.unpack_from("<i", u, 0)[0]
    x_off, y_off = struct.unpack_from("<ii", u, 4)
    use_short = u[12]
    typ = u[13]
    if typ == 1:
        row_count = struct.unpack_from("<h", u, 14)[0]
        p = 16
        xs, ys, vs = [], [], []
        cdt = np.dtype([("x", "<i2"), ("c", "<i2")]) if use_short == 0 else \
            np.dtype([("x", "<i2"), ("c", "<f4")])
        for _ in range(row_count):
            yy, cc = struct.unpack_from("<hh", u, p)
            p += 4
            rec = np.frombuffer(u, cdt, count=cc, offset=p)
            p += cc * cdt.itemsize
            xs.append(rec["x"].astype(np.int64) + x_off)
            ys.append(np.full(cc, yy + y_off, np.int64))
            vs.append(rec["c"].astype(np.float64))
        x = np.concatenate(xs) if xs else np.zeros(0, np.int64)
        y = np.concatenate(ys) if ys else np.zeros(0, np.int64)
        v = np.concatenate(vs) if vs else np.zeros(0)
    elif typ == 2:
        n_pts = struct.unpack_from("<i", u, 14)[0]
        w = struct.unpack_from("<h", u, 18)[0]
        if use_short == 0:
            c = np.frombuffer(u, "<i2", count=n_pts, offset=20)
            keep = c != -32768
        else:
            c = np.frombuffer(u, "<f4", count=n_pts, offset=20)
            keep = ~np.isnan(c)
        idx = np.nonzero(keep)[0]
        row, col = np.divmod(idx, w)
        x = x_off + col.astype(np.int64)
        y = y_off + row.astype(np.int64)
        v = c[idx].astype(np.float64)
    else:
        raise ValueError(f"unknown block type {typ}")
    if len(v) != n_rec and typ == 1:
        raise ValueError(f"record count mismatch {len(v)} vs {n_rec}")
    return x, y, v


# ----------------------------------------------------------------------------
# commands
# ----------------------------------------------------------------------------
def cmd_probe(args):
    h = HicFile(HIC_URL)
    print("version", h.version, "genome", h.genome, "master", h.master)
    print("chroms", h.chroms)
    print("bp res", h.bp_res, "frag res", h.frag_res)
    norms = sorted({(k[0], k[3]) for k in h.norm_index})
    print("norm types x res:", norms)
    print("attributes:", list(h.attributes)[:10])
    print("fetched", h.f.bytes_fetched, "bytes in", h.f.n_requests, "requests")


def cmd_fetch(args):
    os.makedirs(args.cache, exist_ok=True)
    h = HicFile(HIC_URL)
    for name in args.chroms:
        out = os.path.join(args.cache, f"chr{name}_{args.res}.npz")
        if os.path.exists(out):
            print("exists", out)
            continue
        ci = h.chrom_index(name)
        t0 = time.time()
        x, y, v = h.intra_matrix(ci, args.res)
        kr = h.norm_vector("KR", ci, args.res)
        length = h.chroms[ci][1]
        np.savez_compressed(out, x=x.astype(np.int32), y=y.astype(np.int32),
                            v=v.astype(np.float32), kr=kr, length=length,
                            res=args.res, chrom=name)
        print(f"chr{name}: {len(v)} records, sum {v.sum():.0f}, "
              f"{len(kr)} bins, {time.time() - t0:.0f}s, "
              f"total fetched {h.f.bytes_fetched / 1e6:.1f} MB")


def load_oe(path):
    d = np.load(path, allow_pickle=False)
    x, y, v, kr = d["x"], d["y"], d["v"].astype(np.float64), d["kr"]
    nb = len(kr)
    raw = np.zeros((nb, nb))
    np.add.at(raw, (x, y), v)
    # stored as upper triangle (x <= y or x >= y); symmetrise
    off = x != y
    np.add.at(raw, (y[off], x[off]), v[off])
    good = np.isfinite(kr) & (kr > 0)
    rowsum = raw.sum(1)
    good &= rowsum > 0
    keep = np.nonzero(good)[0]
    sub = raw[np.ix_(keep, keep)] / np.outer(kr[keep], kr[keep])
    # integrity check: KR-balanced rows of a complete matrix have equal sums
    rs = sub.sum(1)
    kr_rowsum_cv = float(rs.std() / rs.mean())
    if kr_rowsum_cv > 0.01:
        raise ValueError(f"KR row sums not constant (cv={kr_rowsum_cv:.3g}): "
                         "matrix incomplete or mis-decoded")
    # expected per genomic distance over retained pairs
    n = len(keep)
    dist = np.abs(keep[:, None] - keep[None, :])
    maxd = int(dist.max())
    s = np.bincount(dist.ravel(), weights=sub.ravel(), minlength=maxd + 1)
    c = np.bincount(dist.ravel(), minlength=maxd + 1)
    exp = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        oe = sub / exp[dist]
    oe[~np.isfinite(oe)] = 0.0
    np.fill_diagonal(oe, 0.0)
    return dict(keep=keep, oe=oe, obs=sub, dist=dist, nbins=nb,
                kr_rowsum_cv=kr_rowsum_cv,
                res=int(d["res"]), chrom=str(d["chrom"]),
                length=int(d["length"]), n_raw_records=len(v),
                raw_sum=float(v.sum()))


def knn_edges(score, dist, k):
    """Union of each node's k highest-score partners (score > 0)."""
    n = score.shape[0]
    edges = set()
    for i in range(n):
        s = score[i]
        cand = np.nonzero(s > 0)[0]
        cand = cand[cand != i]
        if len(cand) == 0:
            continue
        # sort by -score, then distance, then index (lexsort: last key primary)
        order = np.lexsort((cand, dist[i, cand], -s[cand]))
        for j in cand[order[:k]]:
            a, b = (i, int(j)) if i < j else (int(j), i)
            edges.add((a, b))
    return sorted(edges)


def read_domains(path, chrom):
    out = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        assert header.startswith("chr1\tx1\tx2"), header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] != chrom:
                continue
            assert p[3] == chrom and p[1] == p[4] and p[2] == p[5]
            out.append((int(p[1]), int(p[2])))
    return sorted(out)


def cmd_build(args):
    import networkx as nx
    summary = []
    for name in args.chroms:
        path = os.path.join(args.cache, f"chr{name}_{args.res}.npz")
        D = load_oe(path)
        keep, res = D["keep"], D["res"]
        pos_of_bin = {int(b): i for i, b in enumerate(keep)}
        # ground-truth groups from Arrowhead domains: bin b is in domain
        # [x1, x2) if its midpoint b*res + res/2 lies in [x1, x2)
        doms = read_domains(args.domains, name)
        groups, seen, n_dup, n_small = [], set(), 0, 0
        for x1, x2 in doms:
            b_lo = int(np.ceil((x1 - res / 2) / res))
            b_hi = int(np.ceil((x2 - res / 2) / res))  # exclusive
            members = tuple(pos_of_bin[b] for b in range(b_lo, b_hi)
                            if b in pos_of_bin)
            if len(members) < 3:
                n_small += 1
                continue
            if members in seen:
                n_dup += 1
                continue
            seen.add(members)
            groups.append(members)
        sizes = [len(g) for g in groups]
        covered = set().union(*map(set, groups)) if groups else set()
        for k in args.ks:
            tag = f"gm12878_chr{name}_{res // 1000}kb_k{k}"
            od = os.path.join(args.out, tag)
            os.makedirs(od, exist_ok=True)
            edges = knn_edges(D["obs"], D["dist"], k)
            n = len(keep)
            with open(os.path.join(od, "edges.tsv"), "w", newline="\n") as fh:
                for a, b in edges:
                    fh.write(f"{a}\t{b}\n")
            with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                for g in groups:
                    fh.write(" ".join(map(str, g)) + "\n")
            with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                fh.write("id\tchrom\tstart\n")
                for i, b in enumerate(keep):
                    fh.write(f"{i}\tchr{name}\t{int(b) * res}\n")
            G = nx.Graph()
            G.add_nodes_from(range(n))
            G.add_edges_from(edges)
            gd = np.array([abs(int(keep[a]) - int(keep[b])) for a, b in edges])
            # edges inside at least one group (diagnostic only)
            memb = [set() for _ in range(n)]
            for gi, g in enumerate(groups):
                for u in g:
                    memb[u].add(gi)
            def within(es):
                return float(np.mean([bool(memb[a] & memb[b]) for a, b in es]))
            intra = within(edges)
            band = knn_edges(1.0 / (1.0 + D["dist"]), D["dist"], k)
            oe_e = knn_edges(D["oe"], D["dist"], k)
            oe_d = np.array([abs(int(keep[a]) - int(keep[b])) for a, b in oe_e])
            row = dict(graph=tag, chrom=f"chr{name}", res=res, k=k, n=n,
                       n_bins_total=D["nbins"], m=len(edges),
                       mean_degree=2 * len(edges) / n,
                       clustering=nx.average_clustering(G),
                       frac_dist_le5=float(np.mean(gd <= 5)),
                       median_dist_bins=float(np.median(gd)),
                       n_components=nx.number_connected_components(G),
                       n_domains_chrom=len(doms), n_groups=len(groups),
                       n_groups_dropped_small=n_small,
                       n_groups_dropped_dup=n_dup,
                       group_size_min=min(sizes), group_size_max=max(sizes),
                       group_size_median=float(np.median(sizes)),
                       frac_nodes_in_group=len(covered) / n,
                       frac_edges_within_group=intra,
                       band_null_m=len(band),
                       band_null_frac_within_group=within(band),
                       oe_rule_median_dist_bins=float(np.median(oe_d)),
                       oe_rule_frac_dist_le5=float(np.mean(oe_d <= 5)),
                       oe_rule_frac_within_group=within(oe_e),
                       kr_rowsum_cv=D["kr_rowsum_cv"],
                       raw_contacts_intra=D["raw_sum"])
            summary.append(row)
            print(json.dumps(row))
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("probe")
    a = sub.add_parser("fetch")
    a.add_argument("--cache", required=True)
    a.add_argument("--chroms", nargs="+", default=["17", "18", "19", "20", "21", "22"])
    a.add_argument("--res", type=int, default=25000)
    a = sub.add_parser("build")
    a.add_argument("--cache", required=True)
    a.add_argument("--domains", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--chroms", nargs="+", default=["17", "18", "19", "20", "21", "22"])
    a.add_argument("--res", type=int, default=25000)
    a.add_argument("--ks", nargs="+", type=int, default=[8, 16, 32])
    args = ap.parse_args()
    {"probe": cmd_probe, "fetch": cmd_fetch, "build": cmd_build}[args.cmd](args)


if __name__ == "__main__":
    main()
