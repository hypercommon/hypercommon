"""score_big.score gives the metrics of cover_metrics on random overlapping covers,
restricted to the truth's nodes as in real_detection_v3.restrict."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cover_metrics as M          # noqa: E402
import score_big as SB             # noqa: E402


def random_cover(rng, ids, k, lo, hi):
    return [set(rng.choice(ids, size=rng.integers(lo, hi), replace=False).tolist()) for _ in range(k)]


def cases():
    rng = np.random.default_rng(0)
    for s in range(40):
        ids = rng.choice(10_000, size=rng.integers(30, 200), replace=False)
        truth = random_cover(rng, ids, rng.integers(2, 12), 3, 25)
        if s % 4 == 0:
            truth += truth[:2]                        # repeated node sets, as SNAP's top 5000
        uni = set().union(*truth)
        outside = rng.choice(np.arange(10_000, 10_050), size=10, replace=False)
        pool = np.r_[ids, outside]
        kind = s % 3
        if kind == 0:                                 # a partition with a giant block
            perm = rng.permutation(pool)
            cut = np.sort(rng.choice(len(perm), size=rng.integers(1, 6), replace=False))
            det = [set(b.tolist()) for b in np.split(perm, cut) if len(b)]
        elif kind == 1:                               # overlapping
            det = random_cover(rng, pool, rng.integers(1, 15), 1, 40)
        else:                                         # overlapping, around one big community
            det = random_cover(rng, pool, rng.integers(1, 8), 1, 15) + [set(pool[: len(pool) * 2 // 3].tolist())]
        yield truth, uni, det


def as_arrays(det):
    det = [sorted(c) for c in det]
    return np.array([u for c in det for u in c], np.int64), np.array([len(c) for c in det], np.int64)


def test_metrics():
    for truth, uni, det in cases():
        T = SB.Truth(None, comms=[np.array(sorted(c)) for c in truth])
        got = SB.score(T, *as_arrays(det))
        cov = [c & uni for c in det]
        cov = [c for c in cov if c]
        cov3 = [c for c in cov if len(c) >= 3]
        dedup = [set(s) for s in {frozenset(c) for c in truth}]
        nodes = sorted(uni)
        rec = M.truth_recovery(cov, truth, nodes)
        sizes = np.array([len(c) for c in truth])
        want = dict(f1=M.f1_best_match(cov, truth, nodes), f1_dedup=M.f1_best_match(cov, dedup, nodes),
                    omega=M.omega(cov, truth, nodes), omega_dedup=M.omega(cov, dedup, nodes),
                    omega_ge3=M.omega(cov3, truth, nodes), coverage=M.coverage(cov, truth),
                    coverage_ge3=M.coverage(cov3, truth), f1_ge3=M.f1_best_match(cov3, truth, nodes),
                    n_comms=len(cov), n_comms_ge3=len(cov3), rec_mean=rec.mean(),
                    rec_small=rec[sizes <= np.median(sizes)].mean())
        for k, v in want.items():
            assert abs(got[k] - v) < 1e-9, (k, got[k], v)


def test_onmi():
    from real_detection_v3 import onmi_large
    for i, (truth, uni, det) in enumerate(cases()):
        if i % 5:
            continue
        T = SB.Truth(None, comms=[np.array(sorted(c)) for c in truth])
        cov = [c & uni for c in det]
        cov = [c for c in cov if c]
        want = onmi_large(cov, truth, sorted(uni)) if cov else 0.0
        assert abs(SB.onmi_of(T, *as_arrays(det)) - want) < 1e-9
