"""
Turn benchmark rows into the two protocols and their uncertainty.

frozen   For each (family, method) the parameter maximising mean Omega over the
         CALIBRATION graphs (all p pooled -- the noise level is not known in
         practice) is frozen, then scored on the TEST graphs.
oracle   For each test graph, the best parameter of each method by Omega. An
         upper bound, labelled as such; same search budget for every method.

Summary per (family, method, protocol):
  curve      mean Omega and ONMI vs p, with 95% bootstrap CIs over test seeds
  area       integrated Omega over p in [0, 0.8] (trapezoid), bootstrap CI;
             paired differences against HyperCommon on the same seeds
  p_cross    first p at which the seed-mean curve drops below 0.5 ("censored"
             if it never does within the grid)
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
B = 2000
rng = np.random.default_rng(12345)


def load(exp):
    df = pd.read_csv(os.path.join(HERE, "results", f"{exp}.csv"))
    df = df[df["error"].isna()].copy()
    # grid extensions ("<method>_x", see methods.GRIDS_X) belong to their base method
    df["method"] = df["method"].str.replace(r"_x$", "", regex=True)
    # SLPA run with 100 iterations is SLPA with one more parameter value
    s100 = df["method"] == "slpa100"
    df.loc[s100, "param"] = df.loc[s100, "param"].astype(str) + "|t100"
    df.loc[s100, "method"] = "slpa"
    df["param"] = df["param"].astype(str)
    df["p"] = df["p"].round(4)
    return df


def frozen_params(df):
    cal = df[df["split"] == "cal"]
    m = cal.groupby(["family", "method", "param"])["omega"].mean().reset_index()
    best = m.loc[m.groupby(["family", "method"])["omega"].idxmax()]
    return {(r.family, r.method): r.param for r in best.itertuples()}


def protocol_table(df, fp):
    test = df[df["split"] == "test"]
    frozen = test[[fp.get((f, m)) == p for f, m, p in zip(test.family, test.method, test.param)]].copy()
    frozen["protocol"] = "frozen"
    idx = test.groupby(["family", "method", "p", "seed"])["omega"].idxmax()
    oracle = test.loc[idx].copy()
    oracle["protocol"] = "oracle"
    return pd.concat([frozen, oracle])


def _area(ps, ys):
    return float(np.trapezoid(ys, ps))


def summarise(tab):
    out = []
    for (fam, meth, prot), g in tab.groupby(["family", "method", "protocol"]):
        piv = g.pivot_table(index="seed", columns="p", values="omega")
        piv_n = g.pivot_table(index="seed", columns="p", values="onmi")
        ps = piv.columns.values
        Y = piv.values                                   # seeds x p
        boot = rng.integers(0, len(Y), size=(B, len(Y)))
        mean = Y.mean(0)
        lo, hi = np.percentile(Y[boot].mean(1), [2.5, 97.5], axis=0)
        areas = np.array([_area(ps, y) for y in Y])
        a_boot = areas[boot].mean(1)
        cross = next((float(p) for p, y in zip(ps, mean) if y < 0.5), None)
        out.append(dict(family=fam, method=meth, protocol=prot,
                        p=ps.tolist(), omega=mean.tolist(), lo=lo.tolist(), hi=hi.tolist(),
                        onmi=piv_n.values.mean(0).tolist(),
                        area=float(areas.mean()), area_lo=float(np.percentile(a_boot, 2.5)),
                        area_hi=float(np.percentile(a_boot, 97.5)),
                        p_cross=cross, seeds=piv.index.tolist(), area_per_seed=areas.tolist()))
    return out


def paired_vs(summ, ref="hc"):
    """Paired bootstrap of area difference (method - ref) on shared seeds."""
    by = {(s["family"], s["method"], s["protocol"]): s for s in summ}
    res = []
    for (fam, meth, prot), s in by.items():
        r = by.get((fam, ref, prot))
        if r is None or meth == ref:
            continue
        common = sorted(set(s["seeds"]) & set(r["seeds"]))
        a = np.array([s["area_per_seed"][s["seeds"].index(x)] for x in common])
        b = np.array([r["area_per_seed"][r["seeds"].index(x)] for x in common])
        d = a - b
        boot = rng.integers(0, len(d), size=(B, len(d)))
        db = d[boot].mean(1)
        res.append(dict(family=fam, method=meth, protocol=prot, diff=float(d.mean()),
                        lo=float(np.percentile(db, 2.5)), hi=float(np.percentile(db, 97.5))))
    return res


def main(exp):
    df = load(exp)
    if exp.startswith("rings"):
        # one noise range for every ring experiment, p in [0, 0.8], applied BEFORE the
        # frozen-parameter calibration (rings_pref also sampled p up to 1.0; those rows
        # stay in the raw file only)
        df = df[df["p"] <= 0.8 + 1e-9]
    fp = frozen_params(df)
    tab = protocol_table(df, fp)
    summ = summarise(tab)
    paired = paired_vs(summ)
    out = {"frozen_params": {f"{k[0]}|{k[1]}": v for k, v in fp.items()},
           "summary": summ, "paired_vs_hc": paired}
    with open(os.path.join(HERE, "results", f"{exp}_summary.json"), "w") as fh:
        json.dump(out, fh)
    # compact console table: area (integrated Omega) per family, both protocols
    rows = [(s["family"], s["method"], s["protocol"], s["area"], s["area_lo"], s["area_hi"], s["p_cross"]) for s in summ]
    t = pd.DataFrame(rows, columns=["family", "method", "protocol", "area", "lo", "hi", "p_cross"])
    for prot in ("frozen", "oracle"):
        print(f"\n=== integrated Omega over p in [0,0.8], protocol={prot} ===")
        w = t[t.protocol == prot].pivot(index="method", columns="family", values="area")
        print(w.round(3).sort_values(w.columns[0], ascending=False).to_string())
    print("\nfrozen params:", {k: v for k, v in fp.items()})


if __name__ == "__main__":
    main(sys.argv[1])
