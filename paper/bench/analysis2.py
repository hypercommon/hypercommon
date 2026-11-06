"""
Analysis of the second benchmark round (results2/<experiment>.csv, run2.py).

For every graph family and method, three ways of choosing the method's parameter:
  oracle      on each test graph, the parameter with the largest Omega
  frozen      one parameter per family, the one with the largest mean Omega over the
              calibration graphs (seeds 1000-1002, all noise levels pooled)
  label_free  on each test graph, the parameter whose cover has the largest
              overlapping modularity (eq) -- the same rule for every method
and, under each, the mean over the noise range (trapezoid over p, divided by its span)
of every metric of the chosen cover, with a 95% bootstrap interval over test seeds,
and the paired difference from HyperCommon under the same protocol.

    python analysis2.py <experiment> [<experiment> ...]
Writes results2/<experiment>_summary2.csv (one row per family x method x protocol x
metric) and results2/<experiment>_paired2.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results2")
METRICS = ["omega", "onmi", "f1", "coverage", "rec_mean", "rec_small", "rec_large", "n_comms"]
B = 2000
rng = np.random.default_rng(12345)


NATIVE_OF = {"leiden_g1": "leiden", "leiden_cpm_dens": "leiden_cpm", "linkcomm_pd": "linkcomm"}
NATIVE_PARAM = {"infomap": 1.0, "slpa": 0.1, "sbm": 0.0}   # reference defaults; sbm: finest level


def select(df):
    """Rows chosen under the protocols, tagged by `protocol`:
    oracle, label_free and frozen as in the module docstring, and
    native      each method's own default or selection rule where it has one (Leiden-RB at
                gamma = 1, Leiden-CPM at gamma = edge density, link communities at maximum
                partition density, Infomap at Markov time 1, SLPA at r = 0.1), and the
                label-free choice otherwise; a selector that is not modularity for partitions."""
    nat_rows = df[df.method.isin(NATIVE_OF)]
    df = df[~df.method.isin(NATIVE_OF)]
    test = df[df.split == "test"]
    key = ["family", "method", "p", "seed"]
    orc = test.loc[test.groupby(key).omega.idxmax()].assign(protocol="oracle")
    lf = test.loc[test.groupby(key).eq.idxmax()].assign(protocol="label_free")
    cal = df[df.split == "cal"].groupby(["family", "method", "param"]).omega.mean().reset_index()
    best = cal.loc[cal.groupby(["family", "method"]).omega.idxmax()][["family", "method", "param"]]
    fro = test.merge(best, on=["family", "method", "param"]).assign(protocol="frozen")
    parts = [orc, lf, fro]
    if True:
        own = nat_rows[nat_rows.split == "test"].assign(method=lambda d: d.method.map(NATIVE_OF))
        fixed = [test[(test.method == m) & (pd.to_numeric(test.param, errors="coerce").round(6) == v)]
                 for m, v in NATIVE_PARAM.items()]
        replaced = set(NATIVE_OF.values()) | set(NATIVE_PARAM)
        have = set(own.method)
        nat = pd.concat([lf[~lf.method.isin(replaced & (have | set(NATIVE_PARAM)))], own] + fixed,
                        ignore_index=True).assign(protocol="native")
        parts.append(nat)
    return pd.concat(parts, ignore_index=True), best


def areas(g, metric):
    """Per test seed, the noise-averaged metric (trapezoid / span)."""
    piv = g.pivot_table(index="seed", columns="p", values=metric)
    if piv.isna().all().all():
        return None
    piv = piv.dropna(axis=1, how="all")
    ps = piv.columns.values.astype(float)
    Y = piv.values
    if len(ps) == 1:
        return pd.Series(Y[:, 0], index=piv.index)
    return pd.Series(np.trapezoid(Y, ps, axis=1) / (ps[-1] - ps[0]), index=piv.index)


def boot(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    b = x[rng.integers(0, len(x), (B, len(x)))].mean(1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)


def main(exp):
    df = pd.read_csv(os.path.join(RES, f"{exp}.csv"), dtype={"param": str}, on_bad_lines="skip")
    extra = os.path.join(RES, f"{exp}_sbm.csv")             # run_sbm.py (graph-tool, under Linux)
    if os.path.exists(extra):
        df = pd.concat([df, pd.read_csv(extra, dtype={"param": str}, on_bad_lines="skip")], ignore_index=True)
    df = df[df.error.isna()].copy()
    df["param"] = df["param"].fillna("")
    sel, frozen_params = select(df)
    summ, paired = [], []
    for (fam, prot), g in sel.groupby(["family", "protocol"]):
        per = {}
        for m, gm in g.groupby("method"):
            for metric in METRICS:
                a = areas(gm, metric)
                if a is None:
                    continue
                per[(m, metric)] = a
                mean, lo, hi = boot(a.values)
                summ.append(dict(experiment=exp, family=fam, method=m, protocol=prot, metric=metric,
                                 mean=mean, lo=lo, hi=hi))
        for (m, metric), a in per.items():
            ref = per.get(("hc", metric))
            if m == "hc" or ref is None:
                continue
            common = a.index.intersection(ref.index)
            d = (a[common] - ref[common]).values
            mean, lo, hi = boot(d)
            paired.append(dict(experiment=exp, family=fam, method=m, protocol=prot, metric=metric,
                               diff=mean, lo=lo, hi=hi))
    pd.DataFrame(summ).to_csv(os.path.join(RES, f"{exp}_summary2.csv"), index=False)
    pd.DataFrame(paired).to_csv(os.path.join(RES, f"{exp}_paired2.csv"), index=False)
    frozen_params.to_csv(os.path.join(RES, f"{exp}_frozen_params2.csv"), index=False)
    sel[["family", "method", "protocol", "p", "seed", "param", "omega", "onmi", "coverage", "n_comms",
         "rec_mean"]].to_csv(os.path.join(RES, f"{exp}_selected2.csv"), index=False)
    print(f"{exp}: {len(summ)} summary rows, {len(paired)} paired rows")


if __name__ == "__main__":
    for e in sys.argv[1:]:
        main(e)
