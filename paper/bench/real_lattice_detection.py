"""
Detection on real networks whose groups are expected to be lattice-like: Hi-C contact
maps with their contact domains (analysis/data_hic/) and protein C-alpha contact
networks with their CATH domains (analysis/data_protein/). Every graph folder holds
edges.tsv (u, v) and groups.txt (one group per line). Methods, grids, seeds and scoring
are those of real_more_detection_v3.py (its task() is reused with this loader).

    python real_lattice_detection.py <hic|protein> [workers] [methods|all] [name substrings]
Output: results/real_lattice_<set>.csv (appended, resumable by (network, method)).
"""
from __future__ import annotations

import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "analysis"))
class _Data(dict):                              # any set <kind> lives in analysis/data_<kind>/
    def __missing__(self, kind):
        return os.path.join(HERE, "..", "analysis", f"data_{kind}")


DATA = _Data()


def load_folder(path):
    import networkx as nx
    G = nx.Graph()
    with open(os.path.join(path, "edges.tsv")) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2 and p[0].lstrip("-").isdigit():
                a, b = int(p[0]), int(p[1])
                if a != b:
                    G.add_edge(a, b)
    truth = []
    with open(os.path.join(path, "groups.txt")) as fh:
        for line in fh:
            g = {int(x) for x in line.split()} & set(G.nodes())
            if len(g) >= 3:
                truth.append(g)
    return G, truth


def _init(kind):
    import real_more
    real_more.load = lambda net: load_folder(os.path.join(DATA[kind], net))
    import methods2
    base = methods2.run_method2

    def run_method2(name, G, ts, seed):
        if name != "sbm":
            yield from base(name, G, ts, seed)
            return
        # nested SBM as run_sbm.py: one fit, every level of the hierarchy a parameter value
        import numpy as np
        import graph_tool.all as gt
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        g = gt.Graph(directed=False)
        g.add_vertex(len(nodes))
        g.add_edge_list([(idx[u], idx[v]) for u, v in G.edges()])
        gt.seed_rng(seed + 1)
        np.random.seed(seed + 1)
        state = gt.minimize_nested_blockmodel_dl(g)
        for level in range(len(state.levels)):
            b = state.project_level(level).get_blocks().a
            groups = {}
            for i, lab in enumerate(b):
                groups.setdefault(int(lab), set()).add(nodes[i])
            if len(groups) <= 1 and level > 0:
                break
            yield level, list(groups.values())

    methods2.run_method2 = run_method2


def main():
    kind = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    import real_more_detection_v3 as R
    methods = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] != "all" else [m for m in R.METHODS if m != "oslom"]
    out = os.path.join(HERE, "results", f"real_lattice_{kind}{os.environ.get('RL_SUFFIX', '')}.csv")
    nets = sorted(d for d in os.listdir(DATA[kind]) if os.path.isfile(os.path.join(DATA[kind], d, "edges.tsv")))
    if len(sys.argv) > 4:                       # optional substrings, any of which a name must contain
        nets = [n for n in nets if any(s in n for s in sys.argv[4].split(","))]
    done = set()
    if os.path.exists(out):
        with open(out, newline="") as fh:
            done = {(r["network"], r["method"]) for r in csv.DictReader(fh) if not r["error"]}
    tasks = [(n, m) for n in nets for m in methods if (n, m) not in done]
    order = {"slpa": 0, "demon": 1, "angel": 2, "walktrap": 3, "cpm": 4}
    tasks.sort(key=lambda t: order.get(t[1], 9))
    print(f"{len(tasks)} tasks on {len(nets)} graphs", flush=True)
    new = not os.path.exists(out)
    with open(out, "a", newline="") as fh, Pool(workers, initializer=_init, initargs=(kind,)) as pool:
        w = csv.DictWriter(fh, fieldnames=R.FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        for i, rows in enumerate(pool.imap_unordered(R.task, tasks), 1):
            w.writerows(rows)
            fh.flush()
            print(f"  {i}/{len(tasks)} {rows[0]['network']} {rows[0]['method']} {rows[0].get('error', '')[:120]}",
                  flush=True)


if __name__ == "__main__":
    main()
