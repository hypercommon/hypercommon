import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
from multiprocessing import Pool
import experiments
from hc_fast import TripleScores
from methods import T_GRID
from self_tune import select

EXP = sys.argv[1]

def job(a):
    fam, p, seed = a
    G, _ = EXPERIMENTS_[fam](p, seed)
    ts = TripleScores(G)
    out = []
    for variant in ("hc", "hc_tri"):
        sel, _ = select(G, ts, variant, T_GRID)
        out += [(fam, p, seed, variant, k, v) for k, v in sel.items()]
    return out

E = experiments.EXPERIMENTS[EXP]
EXPERIMENTS_ = {f: (lambda fam: (lambda p, s: E["build"](fam, p, s)))(fam) for f, fam in E["families"].items()}

if __name__ == "__main__":
    tasks = [(f, p, s) for f in E["families"] for p in E["p_grid"] for s in range(10)]
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as pool:
        rows = [r for rs in pool.imap_unordered(job, tasks) for r in rs]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", f"{EXP}_selftune.csv")
    pd.DataFrame(rows, columns=["family", "p", "seed", "method", "criterion", "t"]).to_csv(out, index=False)
    print("done", len(rows))
