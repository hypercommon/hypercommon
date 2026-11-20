"""
Running time of every method's whole parameter sweep on three benchmark graphs, one
thread, one method at a time on an otherwise idle machine (Table S, Sec. Cost).

Graphs (seed 0, as run2.py builds them): uniform rings (s = 100, k = 16, o_n = 0.05,
p = 0.2), two worlds (o_n = 0.05, p = 0.2) and overlapping LFR (o_n = 0.1, mu = 0.3).
Timed: from the networkx graph to the covers at every value of the method's grid
(methods2.GRIDS2; OSLOM: run_oslom.T_OSLOM; SBM: one fit, its hierarchy levels read
off), wall time; for HyperCommon this includes the scoring of all wedges
(hc_fast.TripleScores). Scoring the covers (Omega, overlapping modularity) is not timed.
Each (graph, method) runs in its own process with OMP/MKL/OpenBLAS threads = 1, three
repetitions; the table reports the median.

    python time_methods.py [methods...]          rows to results/time_methods.csv (appended)
OSLOM and the SBM run under Linux (WSL: the gt environment, OSLOM_BIN set), the others
under Windows, on the same laptop.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "results", "time_methods.csv")
GRAPHS = [("rings_uniform", "k16_s100_o05", 0.2), ("rings_hetero", "two_worlds", 0.2),
          ("lfr", "on10", 0.3)]
METHODS = ["hc", "linkcomm", "cpm", "leiden", "leiden_cpm", "infomap", "lpa", "walktrap",
           "slpa", "demon", "angel", "scan", "scan_j"]
REPS = 3
TIMEOUT = 3600
FIELDS = ["experiment", "family", "p", "method", "rep", "n_nodes", "n_edges", "n_params",
          "wall_s", "platform", "error"]


def child(exp_name, fam_id, p, method):
    from experiments import EXPERIMENTS
    exp = EXPERIMENTS[exp_name]
    G, _ = exp["build"](exp["families"][fam_id], p, 0)
    # every library is imported before the clock starts
    if method == "sbm":
        import graph_tool.all as gt
    else:
        import methods2  # noqa: F401
        from methods2 import HC_LIKE, run_method2
        from hc_fast import TripleScores
        if method not in HC_LIKE + ("oslom",):
            from cdlib import algorithms  # noqa: F401
            import igraph, leidenalg, infomap, cpm_fast  # noqa: F401,E401
    t0 = time.perf_counter()
    if method == "sbm":
        nodes = list(G.nodes())
        idx = {u: i for i, u in enumerate(nodes)}
        g = gt.Graph(directed=False)
        g.add_vertex(len(nodes))
        g.add_edge_list([(idx[u], idx[v]) for u, v in G.edges()])
        gt.seed_rng(1)
        state = gt.minimize_nested_blockmodel_dl(g)
        covers = [state.project_level(lv).get_blocks().a.copy() for lv in range(len(state.levels))]
        n_params = len(covers)
    else:
        ts = TripleScores(G) if method in HC_LIKE else None
        n_params = sum(1 for _ in run_method2(method, G, ts, 0))
    wall = time.perf_counter() - t0
    print(json.dumps(dict(n_nodes=G.number_of_nodes(), n_edges=G.number_of_edges(),
                          n_params=n_params, wall_s=round(wall, 4))))


def main(methods):
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               NUMBA_NUM_THREADS="1")
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        for exp, fam, p in GRAPHS:
            for m in methods:
                for rep in range(REPS):
                    row = dict(experiment=exp, family=fam, p=p, method=m, rep=rep,
                               platform=sys.platform, error="")
                    try:
                        r = subprocess.run([sys.executable, "-u", __file__, "--child", exp, fam, str(p), m],
                                           env=env, capture_output=True, text=True, timeout=TIMEOUT)
                        if r.returncode:
                            row["error"] = r.stderr.strip().splitlines()[-1][:300] if r.stderr else "failed"
                        else:
                            row.update(json.loads(r.stdout.strip().splitlines()[-1]))
                    except subprocess.TimeoutExpired:
                        row["error"] = f"timeout {TIMEOUT} s"
                    w.writerow(row)
                    fh.flush()
                    print(exp, fam, m, rep, row.get("wall_s", row["error"]), flush=True)
                    if row["error"]:
                        break


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        child(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
    else:
        main(sys.argv[1:] or METHODS)
