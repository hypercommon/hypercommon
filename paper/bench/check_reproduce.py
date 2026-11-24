"""Rerun a sample of benchmark tasks under the current Python and compare every scored cover with
the rows stored in results2 (computed under Python 3.12). Used to check that a new Python version
gives the same numbers:

    python check_reproduce.py

Infomap differs from the stored rows on any machine but the one that computed them (its random
numbers depend on the platform; README), so it is reported apart."""
import os, sys, platform
import pandas as pd
B = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, B); os.chdir(B)
import run2

TASKS = [("rings_hetero", "two_worlds", 0.3, 0), ("rings_hetero", "gradient10", 0.1, 3),
         ("rings_uniform", "k8_s50_o05", 0.5, 1), ("lfr", "on10", 0.3, 2)]
METHODS = ["hc", "hc_2test", "hc_tri", "hc_econn", "linkcomm", "cpm", "leiden", "leiden_cpm", "infomap",
           "lpa", "walktrap", "slpa", "demon", "angel", "scan", "scan_j", "bigclam", "egosplit"]
COLS = ["omega", "onmi", "f1", "eq", "n_comms"]
print("python", platform.python_version())
out = []
for exp, fam, p, seed in TASKS:
    old = pd.read_csv(f"results2/{exp}.csv", low_memory=False)
    old = old[(old.family == fam) & (old.p.round(4) == round(p, 4)) & (old.seed == seed)]
    for m in METHODS:
        rows = run2.run_task((exp, fam, p, seed, m))
        new = pd.DataFrame(rows)
        if "error" in new and new.error.astype(str).str.len().gt(0).any():
            out.append(dict(exp=exp, fam=fam, method=m, status="ERROR " + str(new.error.iloc[0])[:80])); continue
        o = old[old.method == m]
        if o.empty:
            out.append(dict(exp=exp, fam=fam, method=m, status="no stored rows")); continue
        key = lambda s: s.map(lambda x: "none" if x is None or (isinstance(x, float) and x != x) or str(x) in ("None", "nan", "") else str(x))
        new["param"] = key(new.param); o = o.assign(param=key(o.param))
        mrg = new.merge(o, on="param", suffixes=("_new", "_old"))
        d = max((mrg[c + "_new"].astype(float) - mrg[c + "_old"].astype(float)).abs().max() for c in COLS)
        out.append(dict(exp=exp, fam=fam, method=m, params=len(new), matched=len(mrg), max_diff=d,
                        status=("same" if len(mrg) == len(new) and d < 1e-9 else "platform (Infomap)" if m == "infomap" else "DIFF")))
R = pd.DataFrame(out)
with pd.option_context("display.width", 200, "display.max_rows", 200):
    print(R.to_string(index=False))
print(R.status.value_counts().to_string())
