"""
Label-free threshold selection for HyperCommon: pick t maximising a quality
function of the cover itself, never the ground truth.

  density   Ahn-Bagrow-Lehmann partition density over communities, using the
            edges inside each community:  D = (2/M) sum_c m_c (m_c-(n_c-1)) /
            ((n_c-2)(n_c-1)); edges outside every community add nothing.
  qov       overlapping modularity with membership-normalised weights
            (Shen et al. / Nicosia et al. family): each node's unit membership
            is split equally across its communities,
            Q = (1/2M) sum_c sum_{i,j in C} (A_ij - k_i k_j / 2M) / (O_i O_j).

`select(G, ts, variant, grid)` returns {criterion: t*}.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def _adj(G):
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    r = [idx[u] for u, v in G.edges()] + [idx[v] for u, v in G.edges()]
    c = [idx[v] for u, v in G.edges()] + [idx[u] for u, v in G.edges()]
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(nodes), len(nodes)))
    return A, idx


def partition_density(A, idx, cover):
    M = A.nnz / 2
    if M == 0 or not cover:
        return 0.0
    total = 0.0
    for com in cover:
        ii = np.fromiter((idx[u] for u in com), int)
        n_c = len(ii)
        if n_c < 3:
            continue
        m_c = A[ii][:, ii].nnz / 2
        total += m_c * (m_c - (n_c - 1)) / ((n_c - 2) * (n_c - 1))
    return 2.0 * total / M


def overlapping_modularity(A, idx, cover):
    n = A.shape[0]
    M2 = A.nnz
    if M2 == 0 or not cover:
        return 0.0
    k = np.asarray(A.sum(axis=1)).ravel()
    O = np.zeros(n)
    for com in cover:
        for u in com:
            O[idx[u]] += 1
    q = 0.0
    for com in cover:
        ii = np.fromiter((idx[u] for u in com), int)
        w = 1.0 / O[ii]
        Ac = A[ii][:, ii]
        q += float(w @ (Ac @ w)) - float((w * k[ii]).sum()) ** 2 / M2
    return q / M2


def select(G, ts, variant, grid):
    A, idx = _adj(G)
    best = {"density": (-np.inf, None), "qov": (-np.inf, None)}
    curves = {"density": [], "qov": []}
    for t in grid:
        cov = ts.communities(variant, t)
        for name, f in (("density", partition_density), ("qov", overlapping_modularity)):
            v = f(A, idx, cov)
            curves[name].append(v)
            if v > best[name][0]:
                best[name] = (v, t)
    return {k: v[1] for k, v in best.items()}, curves
