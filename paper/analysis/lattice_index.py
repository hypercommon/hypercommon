"""
How lattice-like are the recorded groups of a network? For every group, wedges whose
three nodes all lie in the group (the centre adjacent to the other two, inside the
group) are sampled, and the lattice index is the share of them whose three pairs all
have closed-neighbourhood Jaccard at least t = 0.1 -- the statistic of Sec. VII D
(89% of the wedges inside small LFR communities pass, 21% inside large ones). Groups
are weighted by their number of internal wedges. Also reported: the mean Jaccard of the
edges inside groups. Reference values from ring lattices and LFR graphs are computed
the same way.

    python lattice_index.py
Writes figures/data/lattice_index.csv.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "bench"))
OUT = os.path.join(HERE, "..", "figures", "data", "lattice_index.csv")
T = 0.1
PER_GROUP = 400
rng = np.random.default_rng(0)


def jac(nb, a, b):
    A, B = nb[a], nb[b]
    i = len(A & B)
    return i / (len(A) + len(B) - i)


def index(nbrs, groups):
    """nbrs: dict node -> set of neighbours; groups: list of node sets."""
    nb = {u: (s | {u}) for u, s in nbrs.items()}
    num = den = 0.0
    ej = []
    for g in groups:
        g = set(g) & nb.keys()
        if len(g) < 3:
            continue
        cen = [(c, [x for x in nbrs[c] if x in g]) for c in g]
        cen = [(c, L) for c, L in cen if len(L) >= 2]
        w = np.array([len(L) * (len(L) - 1) / 2 for _, L in cen], dtype=float)
        W = w.sum()
        if W == 0:
            continue
        k = int(min(PER_GROUP, W))
        pick = rng.choice(len(cen), size=k, p=w / W)
        ok = 0
        for i in pick:
            c, L = cen[i]
            a, b = rng.choice(len(L), size=2, replace=False)
            a, b = L[a], L[b]
            ok += min(jac(nb, c, a), jac(nb, c, b), jac(nb, a, b)) >= T
        num += W * ok / k
        den += W
        es = [(u, v) for u in g for v in nbrs[u] if v in g and u < v]
        if len(es) > PER_GROUP:
            es = [es[i] for i in rng.choice(len(es), PER_GROUP, replace=False)]
        ej += [jac(nb, u, v) for u, v in es]
    return num / den if den else np.nan, float(np.mean(ej)) if ej else np.nan


def nx_nbrs(G):
    return {u: set(G[u]) - {u} for u in G}


def snap(name):
    S = os.path.join(HERE, "snap")
    nbrs = {}
    with gzip.open(os.path.join(S, f"com-{name}.ungraph.txt.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            a, b = map(int, line.split())
            if a != b:
                nbrs.setdefault(a, set()).add(b)
                nbrs.setdefault(b, set()).add(a)
    f = "top5000" if name != "youtube" else "all"
    with gzip.open(os.path.join(S, f"com-{name}.{f}.cmty.txt.gz"), "rt") as fh:
        groups = {frozenset(map(int, l.split())) for l in fh if l.strip()}
    return nbrs, [g for g in groups if len(g) >= 3]


def main():
    from real_more import ALL, load
    rows = []
    for net in [n for n in ALL if n != "cora"]:
        G, truth = load(net)
        li, mj = index(nx_nbrs(G), truth)
        rows.append(dict(network=net, kind="real", lattice_index=li, mean_edge_J=mj))
        print(rows[-1], flush=True)
    for net in ("amazon", "dblp", "youtube"):
        nbrs, groups = snap(net)
        li, mj = index(nbrs, groups)
        rows.append(dict(network=net, kind="real", lattice_index=li, mean_edge_J=mj))
        print(rows[-1], flush=True)
    import graphs
    import experiments as E
    for p in (0.0, 0.1, 0.2, 0.4):
        G, truth = graphs.rings([50] * 40, [8] * 40, 0.05, p, 0)[:2]
        li, mj = index(nx_nbrs(G), truth)
        rows.append(dict(network=f"rings (50,8)x40, p={p}", kind="reference", lattice_index=li, mean_edge_J=mj))
        print(rows[-1], flush=True)
    for q in (0.1, 0.3, 1.0):
        G, truth = graphs.rings([50] * 40, [8] * 40, 0.05, 0.0, 0, internal=q)[:2]
        li, mj = index(nx_nbrs(G), truth)
        rows.append(dict(network=f"rings (50,8)x40, q={q}", kind="reference", lattice_index=li, mean_edge_J=mj))
        print(rows[-1], flush=True)
    fam = E.EXPERIMENTS["lfr"]["families"]
    for fname in list(fam)[:2]:
        for mu in (0.1, 0.3):
            try:
                G, truth = graphs.lfr(fam[fname], mu, 0)[:2]
            except Exception as e:
                print("lfr failed", fname, mu, e)
                continue
            li, mj = index(nx_nbrs(G), truth)
            rows.append(dict(network=f"LFR {fname}, mu={mu}", kind="reference", lattice_index=li, mean_edge_J=mj))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
