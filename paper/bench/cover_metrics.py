"""
Cover-comparison metrics.

omega   adjusted Omega index (Collins & Dent; Murray et al.), same definition as
        metrics/omega.py, vectorised with sparse co-membership counts. The
        degenerate case Omega_e == 1 (both covers give every pair the same count)
        returns 1.0 only if the two covers agree on every pair, else 0.0.
onmi    overlapping NMI, McDaid-Greene-Hurley (arXiv:1110.2515) "max" variant, as
        implemented by cdlib (overlapping_normalized_mutual_information_MGH).
        Nodes missing from the detected cover are added as singletons first
        (they carry no information), so unassigned nodes are penalised.
coverage fraction of ground-truth nodes covered by at least one detected community.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def _membership(cover, index, n):
    rows, cols = [], []
    for c, com in enumerate(cover):
        for u in com:
            i = index.get(u)
            if i is not None:
                rows.append(i)
                cols.append(c)
    return sp.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(n, max(len(cover), 1)))


def _upper_counts(M):
    C = (M @ M.T).tocoo()
    keep = C.row < C.col
    return C.row[keep].astype(np.int64), C.col[keep].astype(np.int64), C.data[keep].astype(np.int64)


def omega(detected, truth, nodes) -> float:
    nodes = list(nodes)
    n = len(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    P = n * (n - 1) // 2
    r1, c1, d1 = _upper_counts(_membership(truth, index, n))
    r2, c2, d2 = _upper_counts(_membership(detected, index, n))
    k1 = r1 * n + c1
    k2 = r2 * n + c2
    keys = np.union1d(k1, k2)
    v1 = np.zeros(len(keys), np.int64)
    v2 = np.zeros(len(keys), np.int64)
    v1[np.searchsorted(keys, k1)] = d1
    v2[np.searchsorted(keys, k2)] = d2
    zeros = P - len(keys)
    agree = int((v1 == v2).sum()) + zeros
    omega_u = agree / P
    h1 = np.bincount(v1, minlength=1).astype(float)
    h2 = np.bincount(v2, minlength=1).astype(float)
    h1[0] += zeros
    h2[0] += zeros
    L = max(len(h1), len(h2))
    h1 = np.pad(h1, (0, L - len(h1))) / P
    h2 = np.pad(h2, (0, L - len(h2))) / P
    omega_e = float((h1 * h2).sum())
    if omega_e >= 1.0 - 1e-15:
        return 1.0 if agree == P else 0.0
    return (omega_u - omega_e) / (1.0 - omega_e)


def coverage(detected, truth) -> float:
    gt_nodes = set().union(*truth) if truth else set()
    det_nodes = set().union(*detected) if detected else set()
    return len(gt_nodes & det_nodes) / len(gt_nodes) if gt_nodes else 0.0


def _h(p):
    p = np.asarray(p, float)
    out = np.zeros_like(p)
    m = p > 0
    out[m] = -p[m] * np.log2(p[m])
    return out


def _cond_entropy(Ma, Mb, sa, sb, n):
    """
    sum over communities X_i of A of min_j H(X_i | Y_j) with McDaid's
    admissibility rule -- the vectorised form of cdlib's
    coverConditionalEntropy(cover, coverRef, normalized=False).
    """
    I = (Ma.T @ Mb).toarray().astype(float)                 # |X_i & Y_j|
    a = (n - sa[:, None] - sb[None, :] + I) / n
    b = (sb[None, :] - I) / n
    c = (sa[:, None] - I) / n
    d = I / n
    ok = _h(a) + _h(d) > _h(b) + _h(c)
    hy = _h(sb / n) + _h(1 - sb / n)
    hx = _h(sa / n) + _h(1 - sa / n)
    joint = _h(a) + _h(b) + _h(c) + _h(d)
    cond = np.where(ok, joint - hy[None, :], hx[:, None])
    return float(cond.min(axis=1).sum())


def onmi(detected, truth, G) -> float:
    """
    McDaid ONMI ("max" normalisation), identical to cdlib's
    overlapping_normalized_mutual_information_MGH; uncovered nodes are added to
    the detected cover as singletons first, so the node universe is V(G).
    """
    nodes = list(G.nodes())
    n = len(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    covered = set().union(*detected) if detected else set()
    det = [set(c) for c in detected] + [{u} for u in nodes if u not in covered]
    gt = [set(c) for c in truth]
    if not det or not gt:
        return 0.0
    Ma = _membership(det, index, n)
    Mb = _membership(gt, index, n)
    sa = np.asarray(Ma.sum(axis=0)).ravel().astype(float)
    sb = np.asarray(Mb.sum(axis=0)).ravel().astype(float)
    hxy = _cond_entropy(Ma, Mb, sa, sb, n)
    hyx = _cond_entropy(Mb, Ma, sb, sa, n)
    hx = float((_h(sa / n) + _h(1 - sa / n)).sum())
    hy = float((_h(sb / n) + _h(1 - sb / n)).sum())
    return 0.5 * (hx - hxy + hy - hyx) / max(hx, hy)


def onmi_reference(detected, truth, G) -> float:
    from cdlib import NodeClustering
    from cdlib.evaluation import overlapping_normalized_mutual_information_MGH as mgh
    covered = set().union(*detected) if detected else set()
    det = [list(c) for c in detected] + [[u] for u in G.nodes() if u not in covered]
    a = NodeClustering(det, graph=None, method_name="d", overlap=True)
    b = NodeClustering([list(c) for c in truth], graph=None, method_name="t", overlap=True)
    return float(mgh(a, b, normalization="max").score)


def f1_best_match(detected, truth, nodes) -> float:
    """
    Best-match F1 (Yang & Leskovec): average over true communities of the best F1
    with any detected one, and over detected communities of the best F1 with any
    true one; the two averages are averaged. 0 if either cover is empty.
    """
    det = [c for c in detected if c]
    gt = [c for c in truth if c]
    if not det or not gt:
        return 0.0
    nodes = list(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    Md, Mg = _membership(det, index, n), _membership(gt, index, n)
    # sparse: pairs of communities that share no node have F1 = 0
    I = (Md.T @ Mg).tocoo()
    sd = np.asarray(Md.sum(axis=0)).ravel().astype(float)
    sg = np.asarray(Mg.sum(axis=0)).ravel().astype(float)
    F = 2 * I.data.astype(float) / (sd[I.row] + sg[I.col])
    best_d, best_g = np.zeros(len(det)), np.zeros(len(gt))
    np.maximum.at(best_d, I.row, F)
    np.maximum.at(best_g, I.col, F)
    return 0.5 * (best_g.mean() + best_d.mean())


def _f1_best_match_dense(detected, truth, nodes) -> float:
    """Reference implementation of f1_best_match (dense contingency table)."""
    det = [c for c in detected if c]
    gt = [c for c in truth if c]
    if not det or not gt:
        return 0.0
    nodes = list(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    Md, Mg = _membership(det, index, n), _membership(gt, index, n)
    I = (Md.T @ Mg).toarray().astype(float)
    sd = np.asarray(Md.sum(axis=0)).ravel()[:, None]
    sg = np.asarray(Mg.sum(axis=0)).ravel()[None, :]
    F = 2 * I / (sd + sg)
    return 0.5 * (F.max(axis=0).mean() + F.max(axis=1).mean())


def truth_recovery(detected, truth, nodes) -> np.ndarray:
    """
    For every true community, the best Jaccard similarity with any detected
    community (0 if nothing is detected): a community-level score in which a
    small community counts as much as a large one (Omega, being pair-based, is
    dominated by the largest communities).
    """
    gt = [c for c in truth if c]
    det = [c for c in detected if c]
    if not gt:
        return np.zeros(0)
    if not det:
        return np.zeros(len(gt))
    nodes = list(nodes)
    index = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    Md, Mg = _membership(det, index, n), _membership(gt, index, n)
    I = (Md.T @ Mg).toarray().astype(float)
    sd = np.asarray(Md.sum(axis=0)).ravel()[:, None]
    sg = np.asarray(Mg.sum(axis=0)).ravel()[None, :]
    return (I / (sd + sg - I)).max(axis=0)
