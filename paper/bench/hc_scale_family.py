"""
The commonality family at scale: ScaleScores (one wedge enumeration, vectorised)
extended so that every function of hc_family.FUNCTIONS is scored from sparse
matrix products instead of per-pair Python loops.

For a pair (u, v):
  common open neighbours  c  = (A @ A)[u, v]
  resource allocation     ra = (A @ D^-1 @ A)[u, v]
  closed intersection     I  = c + 2 if adjacent else c
  closed sizes            su = d_u + 1, sv = d_v + 1,  U = su + sv - I

`communities_f(f, q)` thresholds the triple scores of f at their q-quantile,
exactly as hc_family.FamilyScores does (test_hc_scale_family.py checks equality).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from hc_family import FUNCTIONS
from hc_scale import ScaleScores


class ScaleFamily(ScaleScores):
    def __init__(self, G=None, edges=None, n=None, block=2_000_000):
        super().__init__(G=G, edges=edges, n=n, block=block)
        n = self.n
        A = self.A.astype(np.float64)
        dinv = sp.diags(1.0 / np.maximum(self.deg, 1))
        R = (A @ dinv @ A).tocsr()
        R.setdiag(0)
        R.eliminate_zeros()
        Rr = np.repeat(np.arange(n), np.diff(R.indptr))
        Rc = R.indices.astype(np.int64)
        m = Rr < Rc
        order = np.argsort(Rr[m] * n + Rc[m])
        self._rkeys = (Rr[m] * n + Rc[m])[order]
        self._rvals = R.data[m][order]

        # the unique pairs of the triples, in the order pair_of indexes them
        t = self.triples.astype(np.int64)
        pk = np.concatenate([t[:, 0] * n + t[:, 1], t[:, 0] * n + t[:, 2], t[:, 1] * n + t[:, 2]])
        uniq = np.unique(pk)
        a, b = uniq // n, uniq % n
        c = self._common(a, b).astype(np.float64)
        adj = self._is_edge(a, b)
        I = c + np.where(adj, 2.0, 0.0)
        su, sv = self.deg[a] + 1.0, self.deg[b] + 1.0
        pos = np.minimum(np.searchsorted(self._rkeys, uniq), len(self._rkeys) - 1)
        ra = np.where(self._rkeys[pos] == uniq, self._rvals[pos], 0.0)
        self.P = {"I": I, "U": su + sv - I, "su": su, "sv": sv, "ra": ra}
        self._fcache = {}

    def triple_score(self, fname):
        s = self._fcache.get(fname)
        if s is None:
            s = np.asarray(FUNCTIONS[fname](self.P), float)[self.pair_of].min(axis=1)
            self._fcache[fname] = s
        return s

    def communities_f(self, fname, q):
        s = self.triple_score(fname)
        thr = np.quantile(s, q) if q > 0 else -np.inf
        saved = self.s3
        self.s3 = s
        try:
            return self.communities(thr)
        finally:
            self.s3 = saved
