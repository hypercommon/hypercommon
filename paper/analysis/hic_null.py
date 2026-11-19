"""
Distance-only control for the Hi-C graphs: on the same retained bins, each bin keeps its
k partners nearest along the chromosome (score 1/(1+|i-j|) in bins, ties to the smaller
id), symmetrised by union, as hic_build.py does with the Hi-C contacts; the groups are the
same contact domains. Written next to each Hi-C graph as <name>_null/.

    python hic_null.py
"""
import os
import shutil

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data_hic")


def main():
    for name in sorted(os.listdir(D)):
        src = os.path.join(D, name)
        if name.endswith("_null") or not os.path.isfile(os.path.join(src, "edges.tsv")):
            continue
        k = int(name.rsplit("_k", 1)[1])
        rows = [l.split("\t") for l in open(os.path.join(src, "nodes.tsv")).read().splitlines()]
        rows = [r for r in rows if r[0].isdigit()]
        ids = np.array([int(r[0]) for r in rows])
        pos = np.array([int(r[2]) for r in rows])
        b = pos // 25000                      # bin index along the chromosome
        assert (ids == np.arange(len(ids))).all() and (np.diff(b) > 0).all()
        E = set()
        n = len(b)
        for i in range(n):
            lo, hi = max(0, i - k), min(n, i + k + 1)
            cand = [j for j in range(lo, hi) if j != i]
            cand.sort(key=lambda j: (abs(b[j] - b[i]), j))
            for j in cand[:k]:
                E.add((min(i, j), max(i, j)))
        dst = os.path.join(D, name + "_null")
        os.makedirs(dst, exist_ok=True)
        with open(os.path.join(dst, "edges.tsv"), "w", newline="\n") as fh:
            fh.writelines(f"{u}\t{v}\n" for u, v in sorted(E))
        for f in ("groups.txt", "nodes.tsv"):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
        print(name, n, len(E))


if __name__ == "__main__":
    main()
