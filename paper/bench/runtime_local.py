"""
The local query of runtime.py timed again after its node cache was made lazy (it had
copied every node's neighbourhood before the search, an O(n + m) step that made the
"local" time grow with n). Same graphs, seed nodes and threshold as runtime.py; run on
an otherwise idle machine. Output: results/runtime_local.csv (family, n, rep, hc_local_s,
community size); analysis/fig_runtime.py uses it in place of runtime.csv's hc_local_s.
"""
import csv
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from runtime import lfr_graph, ring_graph, timed                  # noqa: E402
from hypercommon.algorithm import get_node_community             # noqa: E402
from predicates.jaccard import closed_neighborhood_jaccard_predicate  # noqa: E402

OUT = os.path.join(HERE, "results", "runtime_local.csv")


def main():
    cases = [("rings_k8", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)]
    cases += [("rings_k16", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)]
    cases += [("lfr", n) for n in (1000, 2000, 5000, 10000, 20000, 50000, 100000)]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["family", "n", "rep", "hc_local_s", "size"])
        w.writeheader()
        for fam, n in cases:
            for rep in range(3):
                G = ring_graph(n, 8 if fam == "rings_k8" else 16, rep) if fam.startswith("rings") else lfr_graph(n, rep)
                seed_node = random.Random(rep).choice(list(G.nodes()))
                pred = closed_neighborhood_jaccard_predicate(0.12)
                c, t = min((timed(lambda: get_node_community(G, pred, seed_node)) for _ in range(3)),
                           key=lambda x: x[1])   # best of three, as the query is short
                row = dict(family=fam, n=G.number_of_nodes(), rep=rep, hc_local_s=round(t, 5),
                           size=len(c) if c else 0)
                w.writerow(row)
                fh.flush()
                print(row, flush=True)


if __name__ == "__main__":
    main()
