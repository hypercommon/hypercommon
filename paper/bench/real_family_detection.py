"""
HyperCommon with every commonality function on the full SNAP Amazon and DBLP
graphs (same evaluation as real_detection_v2.py, setting "full"): best-match F1
over the quantile threshold grid for each function; ONMI at each function's best.

    python real_family_detection.py <amazon|dblp>
Output: results/real_family_detection.csv
"""
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cover_metrics as M                                  # noqa: E402
import real_detection_v2 as R                              # noqa: E402
from hc_family import FUNCTIONS                            # noqa: E402
from hc_scale_family import ScaleFamily                    # noqa: E402

OUT = os.path.join(HERE, "results", "real_family_detection.csv")
Q = [round(0.02 * i, 2) for i in range(0, 50)]


def main(net):
    G, truth, dedup = R.load(net)
    universe = set().union(*truth)
    nodes = list(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    e = np.array([(idx[a], idx[b]) for a, b in G.edges()], dtype=np.int64)
    t0 = time.time()
    S = ScaleFamily(edges=e, n=len(nodes))
    S.nodes = nodes
    print(f"{net}: built in {time.time() - t0:.0f}s, {len(S.s3)} triples", flush=True)
    new = not os.path.exists(OUT)
    fh = open(OUT, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=["network", "function", "q", "f1", "f1_dedup", "coverage", "n_comms", "onmi"])
    if new:
        w.writeheader()
    done = set()
    if not new:
        with open(OUT, newline="") as rf:
            done = {(r["network"], r["function"]) for r in csv.DictReader(rf)}
    for f in FUNCTIONS:
        if (net, f) in done:
            continue
        rows = []
        for q in Q:
            cov = R.restrict(S.communities_f(f, q), universe)
            rows.append(dict(network=net, function=f, q=q, f1=M.f1_best_match(cov, truth, universe),
                             f1_dedup=M.f1_best_match(cov, dedup, universe), coverage=M.coverage(cov, truth),
                             n_comms=len(cov), onmi=""))
        b = max(rows, key=lambda r: r["f1"])
        b["onmi"] = R.onmi_large(R.restrict(S.communities_f(f, b["q"]), universe), truth, universe)
        w.writerows(rows)
        fh.flush()
        print(f"  {net} {f:15s} best F1 {b['f1']:.3f} at q={b['q']} ONMI {b['onmi']:.3f}", flush=True)
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1])
