"""
The per-branch threshold (hc_tree.py) on the real networks, exploratory: the task and
scoring of real_more_detection_v3.py, rows to results_tree/.

    python real_tree_detection.py more [workers]          the fourteen Netzschleuder networks
    python real_tree_detection.py <kind> [workers]        a lattice-like kind (real_lattice_detection.py)
"""
from __future__ import annotations

import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "analysis"))
OUTDIR = os.path.join(HERE, "results_tree")


def _init(kind):
    if kind != "more":
        import real_lattice_detection as L
        L._init(kind)
    import methods2
    base = methods2.run_method2

    def run_method2(name, G, ts, seed):
        if name == "hc_tree" and ts is None:
            from hc_fast import TripleScores
            ts = TripleScores(G)
        yield from base(name, G, ts, seed)

    methods2.run_method2 = run_method2


def main():
    kind = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    import real_more_detection_v3 as R
    if kind == "more":
        from real_more import ALL
        nets = [n for n in ALL if n != "cora"]
    else:
        import real_lattice_detection as L
        nets = sorted(d for d in os.listdir(L.DATA[kind]) if os.path.isfile(os.path.join(L.DATA[kind], d, "edges.tsv")))
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"real_{kind}.csv")
    done = set()
    if os.path.exists(out):
        with open(out, newline="") as fh:
            done = {r["network"] for r in csv.DictReader(fh) if not r["error"]}
    tasks = [(n, "hc_tree") for n in nets if n not in done]
    print(f"{kind}: {len(tasks)} graphs", flush=True)
    new = not os.path.exists(out)
    with open(out, "a", newline="") as fh, Pool(workers, initializer=_init, initargs=(kind,)) as pool:
        w = csv.DictWriter(fh, fieldnames=R.FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        for rows in pool.imap_unordered(R.task, tasks):
            w.writerows(rows)
            fh.flush()
            if rows and rows[0].get("error"):
                print("  error", rows[0]["network"], rows[0]["error"][:200], flush=True)


if __name__ == "__main__":
    main()
