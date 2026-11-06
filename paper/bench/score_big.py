"""
Scores the covers of big_job.py (and youtube_job.py's HyperCommon covers of the full
YouTube graph) against the SNAP ground truth, with the rules of real_detection_v4.py:
ground truth = the top-5000 communities (Amazon, DBLP) or all groups (YouTube) of at
least three nodes; detected communities restricted to their nodes (the "universe");
every nonempty restricted community kept, or only those of at least three (*_ge3).

Columns per cover: eq (overlapping modularity on the graph detected on, from the job),
f1 / f1_dedup (best-match F1 against the truth / its distinct node sets), omega /
omega_dedup (adjusted Omega index on the universe nodes), coverage, n_comms, the *_ge3
versions, rec_mean / rec_small / rec_large (best Jaccard of each true community,
averaged over all / the sizes up to / above the median), seconds, build_seconds,
n_wedges, threshold (commonality-function covers). onmi (McDaid, as onmi_large) is
computed on each method's rows of best f1, best omega and largest eq.

Omega never lists node pairs: nodes with the same set of communities are merged into
groups, pair counts are taken between groups (row blocks of S @ S.T), and agreements
only over the pairs that share a true community. test_score_big.py checks equality
with cover_metrics.omega and the other metrics of cover_metrics.

    python score_big.py [--dirs results/big_dl/results_big] [--out results/real_big_scores.csv]
    python score_big.py --summary        # oracle / label-free table from the scores
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import json
import os
import sys
import time
import warnings

import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from real_detection_v3 import onmi_large       # noqa: E402

SNAP = os.path.join(HERE, "..", "analysis", "snap")
YT_FULL = os.path.join(HERE, "results", "youtube_full_dl", "results")
OUT = os.path.join(HERE, "results", "real_big_scores.csv")
SUMMARY = os.path.join(HERE, "results", "real_big_summary.csv")
FIELDS = ["network", "setting", "method", "param", "eq", "f1", "f1_dedup", "omega", "omega_dedup", "onmi",
          "coverage", "n_comms", "f1_ge3", "f1_dedup_ge3", "omega_ge3", "coverage_ge3", "n_comms_ge3",
          "rec_mean", "rec_small", "rec_large", "seconds", "build_seconds", "n_wedges", "threshold"]
WORK = 2e7                                       # group pairs per row block


# ---------------------------------------------------------------------------- truth
class Truth:
    def __init__(self, net, comms=None):
        """Truth of a SNAP network, or of the given communities (arrays of ids)."""
        if comms is None:
            cfile = "all" if net == "youtube" else "top5000"
            comms = []
            with gzip.open(os.path.join(SNAP, f"com-{net}.{cfile}.cmty.txt.gz"), "rt") as fh:
                for line in fh:
                    c = np.unique(np.array(line.split(), np.int64))
                    if len(c) >= 3:
                        comms.append(c)
        comms = [np.unique(np.asarray(c, np.int64)) for c in comms]
        self.uni = np.unique(np.concatenate(comms))
        self.U = len(self.uni)
        self.M = self._matrix(comms)
        distinct = {tuple(c.tolist()) for c in comms}
        self.Md = self._matrix([np.array(c, np.int64) for c in sorted(distinct)])
        self.sizes = np.asarray(self.M.sum(axis=0)).ravel()
        med = np.median(self.sizes)
        self.small, self.large = self.sizes <= med, self.sizes > med
        self.omega_prep = {k: OmegaTruth(m) for k, m in (("raw", self.M), ("dedup", self.Md))}

    def _matrix(self, comms):
        rows = np.concatenate([np.searchsorted(self.uni, c) for c in comms])
        cols = np.repeat(np.arange(len(comms)), [len(c) for c in comms])
        return sp.csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(self.U, len(comms)))

    def restrict(self, node, size):
        """Cover given as SNAP ids grouped by community -> universe x communities matrix
        of the nonempty restricted communities."""
        comm = np.repeat(np.arange(len(size)), size)
        pos = np.minimum(np.searchsorted(self.uni, node), self.U - 1)
        keep = self.uni[pos] == node
        r, c = pos[keep], comm[keep]
        _, c = np.unique(c, return_inverse=True)
        k = int(c.max()) + 1 if len(c) else 0
        return sp.csr_matrix((np.ones(len(r), np.int64), (r, c.ravel())), shape=(self.U, max(k, 1))) if k else None


# ---------------------------------------------------------------------------- metrics
def f1_best(Md, Mg):
    if Md is None or Md.shape[1] == 0 or Mg.shape[1] == 0:
        return 0.0
    I = (Md.T @ Mg).tocoo()
    sd = np.asarray(Md.sum(axis=0)).ravel().astype(float)
    sg = np.asarray(Mg.sum(axis=0)).ravel().astype(float)
    F = 2 * I.data / (sd[I.row] + sg[I.col])
    bd, bg = np.zeros(Md.shape[1]), np.zeros(Mg.shape[1])
    np.maximum.at(bd, I.row, F)
    np.maximum.at(bg, I.col, F)
    return 0.5 * (bg.mean() + bd.mean())


def recovery(Md, Mg):
    best = np.zeros(Mg.shape[1])
    if Md is None:
        return best
    I = (Md.T @ Mg).tocoo()
    sd = np.asarray(Md.sum(axis=0)).ravel().astype(float)
    sg = np.asarray(Mg.sum(axis=0)).ravel().astype(float)
    np.maximum.at(best, I.col, I.data / (sd[I.row] + sg[I.col] - I.data))
    return best


def coverage(Md, U):
    return 0.0 if Md is None else float((np.diff(Md.indptr) > 0).sum()) / U


def groups(M):
    """Nodes with equal community sets merged: group of each node, group x community
    matrix S, group sizes w."""
    M = M.tocsr()
    M.sort_indices()
    nnz = np.diff(M.indptr)
    gid = np.empty(M.shape[0], np.int64)
    one = nnz <= 1
    # rows with at most one community: key = that community (or -1)
    first = np.full(M.shape[0], -1, np.int64)
    first[nnz == 1] = M.indices[M.indptr[:-1][nnz == 1]]
    keys = {}
    for i in np.flatnonzero(~one):
        gid[i] = keys.setdefault(tuple(M.indices[M.indptr[i]:M.indptr[i + 1]].tolist()), len(keys))
    multi = len(keys)
    u, inv = np.unique(first[one], return_inverse=True)
    gid[one] = multi + inv.ravel()
    G = multi + len(u)
    rep = np.zeros(G, np.int64)
    rep[gid] = np.arange(M.shape[0])
    w = np.bincount(gid, minlength=G).astype(np.float64)
    return gid, M[rep], w


def _blocks(S):
    """Row blocks of S with about WORK nonzeros of S @ S.T each."""
    per_col = np.asarray(S.sum(axis=0)).ravel()
    cost = np.cumsum(np.asarray(S @ per_col).ravel())
    start, G = 0, S.shape[0]
    while start < G:
        stop = max(int(np.searchsorted(cost, (cost[start - 1] if start else 0) + WORK, side="right")), start + 1)
        yield start, min(stop, G)
        start = stop


def pair_hist(S, w):
    """Number of node pairs sharing exactly k communities, k >= 1 (index 0 unused)."""
    rs = np.asarray(S.sum(axis=1)).ravel().astype(np.int64)
    h = np.bincount(rs, weights=w * (w - 1) / 2)
    ST = S.T.tocsr()
    for a, b in _blocks(S):
        B = (S[a:b] @ ST).tocoo()
        r = B.row + a
        k = B.col > r
        h = _add(h, np.bincount(B.data[k], weights=w[r[k]] * w[B.col[k]]))
    h[0] = 0
    return h


def _add(x, y):
    L = max(len(x), len(y))
    return np.pad(x, (0, L - len(x))) + np.pad(y, (0, L - len(y)))


class OmegaTruth:
    def __init__(self, Mg):
        self.Mg = Mg.tocsr()
        self.gid, self.S, self.w = groups(self.Mg)
        self.h = pair_hist(self.S, self.w)


def omega(Md, T: OmegaTruth):
    """Adjusted Omega index of the detected cover Md against the truth T on all
    universe nodes (cover_metrics.omega)."""
    U = T.Mg.shape[0]
    P = U * (U - 1) / 2
    if Md is None:
        Md = sp.csr_matrix((U, 1), dtype=np.int64)
    gd, Sd, wd = groups(Md)
    h2 = pair_hist(Sd, wd)
    # joint groups: equal truth set and equal detected set
    jkey = T.gid * (len(wd) + 1) + gd
    uj, jinv = np.unique(jkey, return_inverse=True)
    jinv = jinv.ravel()
    wj = np.bincount(jinv).astype(np.float64)
    rep = np.zeros(len(uj), np.int64)
    rep[jinv] = np.arange(U)
    J1, J2 = T.Mg[rep], Md.tocsr()[rep]
    lab = None
    if np.diff(J2.indptr).max(initial=0) <= 1:     # detected sets of size <= 1: a label
        lab = np.full(len(uj), -1, np.int64)
        nz = np.diff(J2.indptr) == 1
        lab[nz] = J2.indices[J2.indptr[:-1][nz]]
    # pairs sharing a true community: within a joint group, then between joint groups
    v1 = np.asarray(J1.sum(axis=1)).ravel()
    v2 = np.asarray(J2.sum(axis=1)).ravel()
    within = wj * (wj - 1) / 2
    n12 = float(within[(v1 > 0) & (v2 > 0)].sum())
    agree = float(within[(v1 > 0) & (v1 == v2)].sum())
    J1T = J1.T.tocsr()
    for a, b in _blocks(J1):
        B = (J1[a:b] @ J1T).tocoo()
        r = B.row + a
        k = B.col > r
        r, c, x = r[k], B.col[k], B.data[k]
        if lab is not None:
            y = ((lab[r] == lab[c]) & (lab[r] >= 0)).astype(np.int64)
        else:
            y = np.asarray(J2[r].multiply(J2[c]).sum(axis=1)).ravel()
        ww = wj[r] * wj[c]
        n12 += float(ww[y > 0].sum())
        agree += float(ww[x == y].sum())
    h1 = T.h.copy()
    N1, N2 = h1.sum(), h2.sum()
    zeros = P - (N1 + N2 - n12)
    h1[0], h2[0] = P - N1, P - N2
    h1, h2 = _add(h1, 0 * h2), _add(h2, 0 * h1)
    ou = (agree + zeros) / P
    oe = float((h1 * h2).sum()) / P ** 2
    if oe >= 1.0 - 1e-15:
        return 1.0 if agree + zeros >= P else 0.0
    return (ou - oe) / (1.0 - oe)


def score(T: Truth, node, size):
    Md = T.restrict(node, size)
    if Md is not None:
        s = np.asarray(Md.sum(axis=0)).ravel()
        Md3 = Md[:, np.flatnonzero(s >= 3)].tocsr() if (s >= 3).any() else None
    else:
        Md3 = None
    rec = recovery(Md, T.M)
    return dict(f1=f1_best(Md, T.M), f1_dedup=f1_best(Md, T.Md),
                omega=omega(Md, T.omega_prep["raw"]), omega_dedup=omega(Md, T.omega_prep["dedup"]),
                coverage=coverage(Md, T.U), n_comms=0 if Md is None else Md.shape[1],
                f1_ge3=f1_best(Md3, T.M), f1_dedup_ge3=f1_best(Md3, T.Md),
                omega_ge3=omega(Md3, T.omega_prep["raw"]), coverage_ge3=coverage(Md3, T.U),
                n_comms_ge3=0 if Md3 is None else Md3.shape[1],
                rec_mean=float(rec.mean()), rec_small=float(rec[T.small].mean()),
                rec_large=float(rec[T.large].mean()) if T.large.any() else "")


def onmi_of(T: Truth, node, size):
    Md = T.restrict(node, size)
    if Md is None:
        return 0.0
    # onmi_large on membership matrices, uncovered nodes as singletons, in row blocks of
    # about 2e7 cells so that memory stays bounded however many communities there are
    from real_detection_v3 import _cond_entropy_chunked
    from cover_metrics import _h
    n = T.U
    free = np.flatnonzero(np.diff(Md.indptr) == 0)
    S = sp.csr_matrix((np.ones(len(free), np.int64), (free, np.arange(len(free)))), shape=(n, len(free)))
    Ma = sp.hstack([Md, S]).tocsr() if len(free) else Md
    Mb = T.M
    sa = np.asarray(Ma.sum(axis=0)).ravel().astype(float)
    sb = np.asarray(Mb.sum(axis=0)).ravel().astype(float)
    hxy = _cond_entropy_chunked(Ma, Mb, sa, sb, n, chunk=max(20, int(2e7 / Mb.shape[1])))
    hyx = _cond_entropy_chunked(Mb, Ma, sb, sa, n, chunk=max(20, int(2e7 / Ma.shape[1])))
    hx = float((_h(sa / n) + _h(1 - sa / n)).sum())
    hy = float((_h(sb / n) + _h(1 - sb / n)).sum())
    return 0.5 * (hx - hxy + hy - hyx) / max(hx, hy)


# ---------------------------------------------------------------------------- covers
def covers(dirs, youtube_full=True):
    """(network, setting, method, param, (npz file, node key, size key), job row) for
    every cover found."""
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "*.csv"))):
            with open(f, newline="") as fh:
                rows = list(csv.DictReader(fh))
            if not rows or "method" not in rows[0]:
                continue
            for r in rows:
                p = r["param"]
                yield r["network"], r["setting"], r["method"], p, (f[:-4] + ".npz", f"p{p}_node", f"p{p}_size"), r
    if youtube_full and os.path.exists(os.path.join(YT_FULL, "meta.json")):
        meta = json.load(open(os.path.join(YT_FULL, "meta.json")))
        with open(os.path.join(YT_FULL, "summary.csv"), newline="") as fh:
            for r in csv.DictReader(fh):
                t = float(r["t"])
                yield "youtube", "full", r["variant"], str(t),                     (os.path.join(YT_FULL, f"covers_{r['variant']}.npz"), f"t{t:.2f}_node", f"t{t:.2f}_size"),                     dict(eq=r["eq"], seconds=r["seconds"], build_seconds=meta["build_seconds"],
                         n_wedges=meta["wedges_kept"])


_T = None


def _init(net):
    global _T
    _T = Truth(net)


def _load(src):
    f, kn, ks = src
    with np.load(f) as z:
        return z[kn], z[ks]


def _score_task(item):
    setting, method, param, src, job = item
    t1 = time.time()
    row = dict(network=None, setting=setting, method=method, param=param, onmi="",
               **{k: job.get(k, "") for k in ("eq", "seconds", "build_seconds", "n_wedges", "threshold")})
    row.update(score(_T, *_load(src)))
    return row, time.time() - t1


def _onmi_task(item):
    k, src = item
    t1 = time.time()
    return k, onmi_of(_T, *_load(src)), time.time() - t1


def run(dirs, out, youtube_full=True, workers=5):
    from multiprocessing import Pool
    done = set()
    if os.path.exists(out):
        with open(out, newline="") as fh:
            done = {(r["network"], r["setting"], r["method"], r["param"]) for r in csv.DictReader(fh)}
    todo = {}
    for net, setting, method, param, src, job in covers(dirs, youtube_full):
        if (net, setting, method, param) not in done:
            todo.setdefault(net, []).append((setting, method, param, src, job))
    new = not os.path.exists(out)
    fh = open(out, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    for net in sorted(todo):
        print(f"{net}: {len(todo[net])} covers, {workers} workers", flush=True)
        with Pool(workers, initializer=_init, initargs=(net,)) as pool:
            for i, (row, secs) in enumerate(pool.imap_unordered(_score_task, todo[net]), 1):
                row["network"] = net
                w.writerow(row)
                fh.flush()
                print(f"  [{i}/{len(todo[net])}] {row['setting']} {row['method']} {row['param']}: "
                      f"f1 {row['f1']:.3f} omega {row['omega_dedup']:.3f} ({secs:.1f} s)", flush=True)
    fh.close()
    add_onmi(dirs, out, youtube_full, workers)


def add_onmi(dirs, out, youtube_full=True, workers=5):
    """ONMI on each method's rows of best f1, best omega_dedup and largest eq."""
    from multiprocessing import Pool
    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh))
    pick = set()
    for col in ("f1", "omega_dedup", "eq"):
        best = {}
        for r in rows:
            k = (r["network"], r["setting"], r["method"])
            if r[col] != "" and (k not in best or float(r[col]) > float(best[k][col])):
                best[k] = r
        pick |= {(r["network"], r["setting"], r["method"], r["param"]) for r in best.values() if r["onmi"] == ""}
    if not pick:
        return
    got = {}
    for net in sorted({p[0] for p in pick}):
        items = [((n2, s, m, p), src) for n2, s, m, p, src, _ in covers(dirs, youtube_full)
                 if n2 == net and (n2, s, m, p) in pick]
        with Pool(workers, initializer=_init, initargs=(net,)) as pool:
            for i, (k, v, secs) in enumerate(pool.imap_unordered(_onmi_task, items), 1):
                got[k] = v
                print(f"  onmi {k}: {v:.3f} ({secs:.0f} s)", flush=True)
                if i % 10 == 0 or i == len(items):          # keep what is done
                    _write_onmi(rows, got, out)


def _write_onmi(rows, got, out):
    for r in rows:
        k = (r["network"], r["setting"], r["method"], r["param"])
        if k in got:
            r["onmi"] = got[k]
    with open(out + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    os.replace(out + ".tmp", out)


# ---------------------------------------------------------------------------- summary
def summary(out, path):
    """Per network, setting and method: the row of best f1, of best omega_dedup (oracle)
    and of largest eq (label-free). 'hc_family' takes the commonality functions together
    (its label-free row: largest eq over every function and quantile)."""
    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        if r["method"].startswith("hcf_"):
            rows.append(dict(r, method="hc_family", param=f"{r['method'][4:]}@{r['param']}"))
    res = []
    keys = sorted({(r["network"], r["setting"], r["method"]) for r in rows})
    for k in keys:
        rs = [r for r in rows if (r["network"], r["setting"], r["method"]) == k]
        for rule, col in (("oracle_f1", "f1"), ("oracle_omega", "omega_dedup"), ("label_free", "eq")):
            ok = [r for r in rs if r[col] != ""]
            if ok:
                res.append(dict(rule=rule, **max(ok, key=lambda r: float(r[col]))))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rule"] + FIELDS)
        w.writeheader()
        w.writerows(res)
    print(f"{len(res)} rows -> {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", default=os.path.join(HERE, "results", "big_dl", "results_big"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--no-youtube-full", action="store_true", help="skip youtube_job.py's covers")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--onmi-only", action="store_true", help="only add the missing ONMI values")
    a = ap.parse_args()
    if a.summary:
        summary(a.out, SUMMARY if a.out == OUT else a.out[:-4] + "_summary.csv")
    elif a.onmi_only:
        add_onmi(a.dirs.split(","), a.out, not a.no_youtube_full, a.workers)
    else:
        run(a.dirs.split(","), a.out, not a.no_youtube_full, a.workers)
