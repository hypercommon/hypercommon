"""
Exact (rational) recovery windows of HyperCommon on small ring-lattice configurations,
used to check the statements of Sec. V D and the remarks after Theorem 2 and
Proposition 4.

A configuration is a list of rings (size, degree) and a list of merges, each a tuple
of (ring, position) items that are identified into one node. The cover is constant on
every interval (v_i, v_{i+1}] between consecutive distinct commonality values, so it
is evaluated once per interval; the result is the union of the intervals on which the
cover equals the ground truth.

    python exact_windows.py        # runs every check and prints PASS/FAIL lines
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations


def build(rings, merges):
    """Closed neighbourhoods and ground truth of the merged ring lattices."""
    ident = {}
    for grp in merges:
        rep = grp[0]
        for item in grp:
            ident[item] = rep
    node = lambda r, i: ident.get((r, i), (r, i))
    nbr = {}
    truth = []
    for r, (s, k) in enumerate(rings):
        comm = set()
        for i in range(s):
            u = node(r, i)
            comm.add(u)
            for d in range(1, k // 2 + 1):
                v = node(r, (i + d) % s)
                nbr.setdefault(u, set()).add(v)
                nbr.setdefault(v, set()).add(u)
        truth.append(frozenset(comm))
    return nbr, set(truth)


def scores(nbr):
    cl = {u: nb | {u} for u, nb in nbr.items()}
    J = {}

    def j(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in J:
            A, B = cl[a], cl[b]
            J[key] = F(len(A & B), len(A | B))
        return J[key]

    tri = {}
    for c, nb in nbr.items():
        for a, b in combinations(sorted(nb), 2):
            key = tuple(sorted((a, b, c)))
            if key not in tri:
                tri[key] = min(j(a, b), j(a, c), j(b, c))
    return tri


def cover(tri, t):
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    adm = [T for T, sc in tri.items() if sc >= t]
    for T in adm:
        a, b, c = T
        p = [(a, b), (a, c), (b, c)]
        for q in p[1:]:
            ra, rb = find(p[0]), find(q)
            if ra != rb:
                parent[ra] = rb
    comps = {}
    for T in adm:
        comps.setdefault(find((T[0], T[1])), set()).update(T)
    return {frozenset(c) for c in comps.values()}


def window(rings, merges):
    """Recovery set as a list of (lo, hi] intervals (lo = None for 'from 0')."""
    nbr, truth = build(rings, merges)
    tri = scores(nbr)
    vals = sorted(set(tri.values()))
    ok, prev = [], None
    for v in vals:
        if cover(tri, v) == truth:
            ok.append((prev, v))
        prev = v
    merged = []
    for lo, hi in ok:
        if merged and merged[-1][1] == lo:
            merged[-1] = (merged[-1][0], hi)
        else:
            merged.append((lo, hi))
    return merged


def show(w):
    return " U ".join(f"({lo}, {hi}]" for lo, hi in w) if w else "empty"


def check(name, got, want):
    print(("PASS " if got == want else "FAIL ") + name + f": {show(got)}" + ("" if got == want else f"  want {show(want)}"))


def two_ring(sa, ka, sb, kb, pos=((0, 0),)):
    return [(sa, ka), (sb, kb)], [((0, a), (1, b)) for a, b in pos]


if __name__ == "__main__":
    # Eq. (5): one isolated merge between rings of degrees ka, kb
    def U(x, y):
        return min(F(x, x + y + 2), max(F(x - 2, x + 4), F(x - 1, x + y + 3)))
    for ka, kb in [(8, 4), (16, 8), (32, 8), (12, 4)]:
        r, m = two_ring(4 * ka + 10, ka, 4 * kb + 10, kb)
        check(f"Eq.5 ({ka},{kb})", window(r, m), [(F(1, ka + kb + 1), min(U(ka, kb), U(kb, ka)))])

    # Mixed families: intersection of the windows of their merge types (three rings,
    # merges a-b of type (ka,kb) and b-c of type (kb,kb), far apart)
    for ka, kb, want in [(32, 8, [(F(1, 17), F(4, 21))]), (16, 8, [(F(1, 17), F(4, 13))]),
                         (28, 4, [(F(1, 9), F(2, 17))]), (30, 4, []), (32, 4, [])]:
        rings = [(4 * ka + 10, ka), (8 * kb + 20, kb), (4 * kb + 10, kb)]
        merges = [((0, 0), (1, 0)), ((1, 4 * kb + 10), (2, 0))]
        check(f"mixed ({ka},{kb})+({kb},{kb})", window(rings, merges), want)

    # Two rings sharing two merged nodes at ring distances DA, DB (k = 8, s = 60):
    # fused (no recovery) for small t; the region and example of Proposition 4
    for DA, DB, want in [(11, 2, [(F(3, 23), F(4, 9))]), (4, 4, None)]:
        r, m = two_ring(60, 8, 60, 8, pos=((0, 0), (DA, DB)))
        w = window(r, m)
        if want is None:
            print(f"INFO two shared nodes {DA},{DB} apart (k=8): {show(w)}")
        else:
            check(f"two shared nodes {DA},{DB} apart (k=8)", w, want)
    # region map for k = 8: is the lower edge raised above 1/17?
    k = 8
    bad = []
    for DA in range(1, 3 * k):
        for DB in range(1, DA + 1):
            r, m = two_ring(70, k, 70, k, pos=((0, 0), (DA, DB)))
            w = window(r, m)
            raised = not w or w[0][0] != F(1, 2 * k + 1)
            region = (DA <= k and DB <= k) or (min(DA, DB) <= k // 2 and max(DA, DB) <= 3 * k // 2)
            if raised != region:
                bad.append((DA, DB, show(w)))
    print(("PASS" if not bad else "FAIL") + f" fusion region k=8 (both<=k, or min<=k/2 and max<=3k/2): {len(bad)} mismatches {bad[:5]}")

    # Three rings pairwise sharing single nodes that form a wedge (shared nodes 2 apart)
    for k, want in [(8, [(F(4, 13), F(4, 9))]), (4, [])]:
        s = 6 * k
        rings = [(s, k)] * 3
        merges = [((0, 0), (1, 0)), ((0, 2), (2, 0)), ((1, 2), (2, 2))]
        check(f"three rings, shared nodes 2 apart (k={k})", window(rings, merges), want)

    # Separation: shared nodes k apart on both rings raise the lower edge to 1/k
    for k in (4, 6, 8):
        r, m = two_ring(6 * k, k, 6 * k, k, pos=((0, 0), (k, k)))
        w = window(r, m)
        print(("PASS" if w and w[0][0] == F(1, k) else "FAIL") + f" separation k={k}: {show(w)}")
        r, m = two_ring(6 * k, k, 6 * k, k, pos=((0, 0), (k + 1, k + 1)))
        w = window(r, m)
        print(("PASS" if w and w[0][0] == F(1, 2 * k + 1) else "FAIL") + f" separation k+1, k={k}: {show(w)}")

    # Small rings, one merge: the window of Eq. (4) down to s = k+4 (k = 4, 6), s = k+2 (k = 8..12)
    def eq4(k):
        J3, Ju1, Ju2 = F(k - 2, k + 4), F(k, 2 * k + 2), F(k - 1, 2 * k + 3)
        return [(F(1, 2 * k + 1), min(Ju1, max(J3, Ju2)))]
    for k in (4, 6, 8, 10, 12):
        smin = k + 4 if k <= 6 else k + 2
        res = []
        for s in range(k + 2, 2 * k + 3):
            w = window(*two_ring(s, k, s, k))
            res.append((s, w == eq4(k)))
        first = min(s for s, ok in res if all(o for s2, o in res if s2 >= s))
        print(("PASS" if first == smin else "FAIL") + f" small rings k={k}: Eq.4 holds for s >= {first} (claim {smin})")
    for s in (6, 7, 9):
        print(f"INFO k=4 s={s}: {show(window(*two_ring(s, 4, s, 4)))}")

    # A node shared by m rings (C6)
    for k, ms in [(4, (2, 3, 4)), (6, (3,)), (8, (3, 4))]:
        for mm in ms:
            rings = [(4 * k + 10, k)] * mm
            merges = [tuple((r, 0) for r in range(mm))]
            w = window(rings, merges)
            up = min(F(k, mm * k + 2), max(F(k - 2, k + 4), F(k - 1, mm * k + 3)))
            want = [(F(1, 2 * k + 1), up)] if F(1, 2 * k + 1) < up else []
            check(f"node shared by m={mm} rings (k={k})", w, want)
