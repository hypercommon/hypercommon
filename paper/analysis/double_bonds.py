"""
Double bonds in the ring benchmarks. Two rings that share two merged nodes u1, u2
lying within ring distance k of each other on BOTH rings have the pair {u1, u2} inside
a wedge of each ring, so the two rings' triples can meet in that pair without any
triple spanning both rings (Proposition "shared pairs" in the paper). For every ring
family at p = 0, test seeds 0-9: the number of ring pairs that share nodes, how many
of them are double-bonded, and the share of rings involved in at least one double bond.
Writes figures/data/double_bonds.csv.
"""
import os
import sys
from itertools import combinations

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench"))
import experiments  # noqa: E402


def analyse(fam):
    sizes, ks = fam["sizes"], fam["ks"]
    offs = [sum(sizes[:r]) for r in range(len(sizes))]
    rows = []
    for seed in range(10):
        import random
        import math
        import graphs
        from generators.ring_lattice import ring_lattice
        from generators.overlap import apply_overlap
        rng = random.Random(f"rings|{sizes}|{ks}|{fam['overlap']}|0.0|{seed}")
        G = ring_lattice(list(sizes), list(ks))
        merged = apply_overlap(G, list(sizes), fam["overlap"], rng)   # {absorbed v: surviving u}
        ring_of = [r for r, s in enumerate(sizes) for _ in range(s)]
        # every merged node: its position on its own ring and on the absorbed node's ring
        pos = {}
        for v, u in merged.items():
            pos.setdefault(u, {})[ring_of[u]] = u - offs[ring_of[u]]
            pos[u][ring_of[v]] = v - offs[ring_of[v]]
        by_pair = {}
        for u, p in pos.items():
            if len(p) == 2:
                by_pair.setdefault(tuple(sorted(p)), []).append(u)
        dbl = set()
        for (ra, rb), us in by_pair.items():
            for x, y in combinations(us, 2):
                da = abs(pos[x][ra] - pos[y][ra]); da = min(da, sizes[ra] - da)
                db = abs(pos[x][rb] - pos[y][rb]); db = min(db, sizes[rb] - db)
                if da <= ks[ra] and db <= ks[rb]:
                    dbl.add((ra, rb))
        involved = {r for pr in dbl for r in pr}
        rows.append(dict(seed=seed, merges=len(merged), ring_pairs_sharing=len(by_pair),
                         double_bonded_pairs=len(dbl), rings=len(sizes),
                         rings_in_double_bond=len(involved), rings_fused=fused_rings(fam, seed)))
    return pd.DataFrame(rows)


def fused_rings(fam, seed):
    """Rings that exact HyperCommon fuses with another ring on the noise-free graph, at a
    threshold just above the lower edge 1/(2 k_min + 1) of the family's recovery window
    (by nesting, every fusion inside the window shows here): a ring is fused when one
    detected community holds at least three of its own unmerged nodes and at least three
    of another ring's."""
    import graphs
    from hc_fast import TripleScores
    G, truth = graphs.rings(fam["sizes"], fam["ks"], fam["overlap"], 0.0, seed)
    shared = {u for i, a in enumerate(truth) for b in truth[i + 1:] for u in a & b}
    own = [c - shared for c in truth]
    t = 1 / (2 * min(fam["ks"]) + 1) + 1e-9
    fused = set()
    for C in TripleScores(G).communities("hc", t):
        hit = [r for r, o in enumerate(own) if len(C & o) >= 3]
        if len(hit) >= 2:
            fused.update(hit)
    return len(fused)


def main():
    out = []
    for exp in ("rings_uniform", "rings_hetero", "rings_hetero_overlap"):
        for name, fam in experiments.EXPERIMENTS[exp]["families"].items():
            if fam["overlap"] == 0:
                continue
            d = analyse(fam)
            out.append(dict(experiment=exp, family=name, **d.mean(numeric_only=True).round(2).to_dict()))
    d = pd.DataFrame(out).drop(columns="seed")
    d["share_rings_double_bonded"] = (d.rings_in_double_bond / d.rings).round(3)
    d["share_rings_fused"] = (d.rings_fused / d.rings).round(3)
    d.to_csv(os.path.join(HERE, "..", "figures", "data", "double_bonds.csv"), index=False)
    print(d.to_string())


if __name__ == "__main__":
    main()
