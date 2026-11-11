"""
OSLOM (OSLOM2 binary, as run_oslom.py) on SNAP Amazon, DBLP and YouTube, full graphs and
the subgraphs induced by the ground-truth nodes: covers in big_job.py's format, scored
afterwards with the other methods by score_big.py. One task per (network, setting,
threshold t of T_OSLOM), seed 1 as big_job.py's Leiden and Infomap; a run that exceeds
the time limit is recorded in DIR/errors.txt and has no cover.

    python run_oslom_big.py [--workers 6] [--only amazon,dblp,youtube] [--settings full,induced]
                            [--hours 48] [--out results/big_dl/results_big]
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def task(args):
    net, setting, t, out, hours = args
    from big_job import Store, adjacency, load, to_arrays
    from run_oslom import oslom_cover_edges
    import numpy as np
    st = Store(out, net, setting, "oslom", f"_{t}")
    if st.done():
        return f"{time.strftime('%H:%M:%S')} skip  {net} {setting} t={t}"
    t0 = time.time()
    workdir = tempfile.mkdtemp(prefix="oslom_")
    try:
        e, ids = load(net, setting)
        A = adjacency(e, len(ids))
        cover, secs = oslom_cover_edges(e.tolist(), t, 1, workdir, timeout=hours * 3600)
        st.add(t, *to_arrays(cover), A, ids, secs)
        st.save()
        return f"{time.strftime('%H:%M:%S')} done  {net} {setting} t={t} ({secs / 60:.1f} min, {len(cover)} groups)"
    except subprocess.TimeoutExpired:
        msg = f"{net} {setting} t={t}: no cover within {hours:g} h"
    except Exception:
        msg = f"{net} {setting} t={t}\n{traceback.format_exc()}"
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    with open(os.path.join(out, "errors.txt"), "a") as fh:
        fh.write(f"oslom {msg}\n")
    return f"{time.strftime('%H:%M:%S')} ERROR {msg.splitlines()[0]} ({(time.time() - t0) / 60:.1f} min)"


def main():
    from run_oslom import T_OSLOM
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", default="amazon,dblp,youtube")
    ap.add_argument("--settings", default="full,induced")
    ap.add_argument("--hours", type=float, default=48.0)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "big_dl", "results_big"))
    a = ap.parse_args()
    tasks = [(net, s, t, a.out, a.hours) for net in a.only.split(",") for s in a.settings.split(",")
             for t in T_OSLOM]
    tasks.sort(key=lambda x: (x[0] != "youtube", x[1] != "full"))     # the longest first
    print(f"{len(tasks)} tasks, {a.workers} workers, limit {a.hours:g} h per run", flush=True)
    with Pool(a.workers, maxtasksperchild=1) as pool:
        for i, m in enumerate(pool.imap_unordered(task, tasks), 1):
            print(f"[{i}/{len(tasks)}] {m}", flush=True)
    print(f"{time.strftime('%H:%M:%S')} ALL DONE", flush=True)


if __name__ == "__main__":
    main()
