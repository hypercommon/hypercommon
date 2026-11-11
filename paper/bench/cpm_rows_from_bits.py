"""
Clique percolation rows of real_more_v2.csv and real_more_v3.csv for a network on which
bench/cpm_fast.py is too slow (facebook_107), from the covers that cpm_bits.py saved:

    python cpm_bits.py facebook_107 3 4 5 6 7 --save covers.pkl      (numba env)
    python cpm_rows_from_bits.py facebook_107 covers.pkl              (hc env)

The two drivers' own task() functions score the covers exactly as in their runs; only
the clique percolation itself is replaced by the saved covers (cpm_bits.py gives the
output of networkx's k_clique_communities, as cpm_fast.py does). In real_more_v3.csv the
first k carries the whole percolation time, as there.
"""
from __future__ import annotations

import csv
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "analysis"))


def main():
    net, path = sys.argv[1], sys.argv[2]
    with open(path, "rb") as f:
        saved = pickle.load(f)
    covers = {int(k): [frozenset(c) for c in v] for k, v in saved["covers"].items()}

    import cpm_fast

    def from_saved(G, ks=(3, 4, 5, 6, 7, 8), cliques=None):
        return {k: covers[k] for k in ks}
    cpm_fast.k_clique_covers = from_saved

    import real_more_detection_v2 as V2
    import real_more_detection_v3 as V3
    for mod in (V2, V3):
        with open(mod.OUT, newline="") as fh:
            if any(r["network"] == net and r["method"] == "cpm" and not r["error"] for r in csv.DictReader(fh)):
                print(f"{os.path.basename(mod.OUT)}: {net} cpm already there")
                continue
        rows = mod.task((net, "cpm"))
        bad = [r for r in rows if r.get("error")]
        if bad:
            raise SystemExit(bad[0]["error"])
        if mod is V3:
            rows[0]["seconds"] = round(saved["secs"], 3)
        with open(mod.OUT, "a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=mod.FIELDS, extrasaction="ignore").writerows(rows)
        print(f"{os.path.basename(mod.OUT)}: {len(rows)} rows, k = {[r['param'] for r in rows]}")


if __name__ == "__main__":
    main()
