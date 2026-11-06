"""
Derived result files that the tables read, produced from the benchmark summaries
(bench/analysis.py output). Run after bench/analysis.py and before make_tables.py.

  figures/data/family_vs_jaccard_paired.csv   each commonality function minus Jaccard,
                                              paired by seed, with a bootstrap 95% CI
  figures/data/hetero_by_overlap.csv          unequal rings at o_n = 0, 0.05, 0.10:
                                              HC, Leiden, best other method, and the
                                              paired comparisons used in the overlap table
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "bench", "results")
DATA = os.path.join(HERE, "..", "figures", "data")


def family_vs_jaccard():
    S = json.load(open(os.path.join(RES, "commonality_family_summary.json")))
    by = {(x["family"], x["method"]): x for x in S["summary"] if x["protocol"] == "oracle"}
    rng = np.random.default_rng(1)
    rows = []
    for fam in sorted({k[0] for k in by}):
        ref = by[(fam, "hcf_jaccard")]
        sp = ref["p"][-1] - ref["p"][0]
        for m in ["hcf_cosine", "hcf_max_norm", "hcf_min_norm", "hcf_i2_over_u", "hcf_count", "hcf_resource_alloc"]:
            x = by[(fam, m)]
            common = sorted(set(x["seeds"]) & set(ref["seeds"]))
            d = np.array([x["area_per_seed"][x["seeds"].index(s)] - ref["area_per_seed"][ref["seeds"].index(s)]
                          for s in common]) / sp
            b = d[rng.integers(0, len(d), (2000, len(d)))].mean(1)
            rows.append(dict(family=fam, method=m[4:], diff=d.mean(), lo=np.percentile(b, 2.5),
                             hi=np.percentile(b, 97.5)))
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "family_vs_jaccard_paired.csv"), index=False)


def hetero_by_overlap():
    rows = []
    for exp in ("rings_hetero", "rings_hetero_overlap"):
        S = json.load(open(os.path.join(RES, f"{exp}_summary.json")))
        span = {x["family"]: x["p"][-1] - x["p"][0] for x in S["summary"]}
        P = pd.DataFrame(S["paired_vs_hc"])
        P = P[P.protocol == "oracle"]
        summ = pd.DataFrame([dict(family=x["family"], method=x["method"], area=x["area"])
                             for x in S["summary"] if x["protocol"] == "oracle"])
        for fam, g in summ.groupby("family"):
            base = fam.replace("_o00", "").replace("_o10", "")
            if base not in ("two_worlds", "three_big_many_small", "gradient10", "sizes_only_k8"):
                continue
            g = g.set_index("method").area / span[fam]
            oth = g.drop([m for m in ("hc", "hc_tri", "hc_2test") if m in g.index])
            bo = oth.idxmax()
            r = P[(P.family == fam) & (P.method == bo)].iloc[0]
            t = P[(P.family == fam) & (P.method == "hc_tri")].iloc[0]
            state = lambda x: "ahead" if x["hi"] < 0 else ("behind" if x["lo"] > 0 else "tie")
            o = 0 if fam.endswith("_o00") else (10 if fam.endswith("_o10") else 5)
            rows.append(dict(family=base, overlap=o, hc=round(g["hc"], 3), best_other=bo,
                             best_other_val=round(oth.max(), 3), vs_best=state(r), tri=round(g["hc_tri"], 3),
                             vs_tri=state(t), leiden=round(g["leiden"], 3)))
    pd.DataFrame(rows).sort_values(["family", "overlap"]).to_csv(os.path.join(DATA, "hetero_by_overlap.csv"),
                                                                    index=False)


def selftune_vs_oracle():
    """Label-free threshold (max overlapping modularity) against the oracle threshold and
    against every other method under the frozen protocol, as mean Omega over the noise range."""
    rows = []
    for exp in ("rings_uniform", "rings_hetero", "rings_pref", "lfr"):
        st = pd.read_csv(os.path.join(RES, f"{exp}_selftune.csv"))
        st = st[(st.method == "hc") & (st.criterion == "qov")]
        raw = pd.read_csv(os.path.join(RES, f"{exp}.csv"), usecols=["family", "p", "seed", "split", "method",
                                                                     "param", "omega", "error"])
        raw = raw[(raw.split == "test") & (raw.method == "hc") & raw.error.isna()].copy()
        if exp.startswith("rings"):
            raw, st = raw[raw.p <= 0.8 + 1e-9], st[st.p <= 0.8 + 1e-9]
        for d in (raw, st):
            d["p"] = d.p.round(3)
        raw["t"] = pd.to_numeric(raw.param, errors="coerce").round(4)
        st = st.assign(t=st.t.round(4))
        oracle = raw.groupby(["family", "p", "seed"]).omega.max().rename("oracle")
        sel = st.merge(raw, on=["family", "p", "seed", "t"]).set_index(["family", "p", "seed"]).omega.rename("self")
        assert len(sel) == len(st), (exp, len(sel), len(st))
        m = pd.concat([oracle, sel], axis=1).dropna().groupby(["family", "p"]).mean().reset_index()
        S = json.load(open(os.path.join(RES, f"{exp}_summary.json")))
        frozen = {(x["family"], x["method"]): x["area"] / (x["p"][-1] - x["p"][0]) for x in S["summary"]
                  if x["protocol"] == "frozen" and x["method"] not in ("hc", "hc_tri", "hc_2test")}
        for fam, g in m.groupby("family"):
            ps = g.p.values
            a_self, a_or = np.trapezoid(g["self"], ps), np.trapezoid(g["oracle"], ps)
            fo = {k[1]: v for k, v in frozen.items() if k[0] == fam}
            bo = max(fo, key=fo.get)
            rows.append(dict(experiment=exp, family=fam, self_tuned=a_self / (ps[-1] - ps[0]),
                             oracle=a_or / (ps[-1] - ps[0]), ratio=a_self / a_or,
                             best_other_frozen=bo, best_other_frozen_omega=fo[bo]))
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "selftune_vs_oracle.csv"), index=False)


def real_label_free():
    """F1 of HyperCommon at its label-free threshold over F1 at its best threshold, on the
    sixteen real networks (full SNAP graphs and the Netzschleuder networks), both scoring rules."""
    more = pd.read_csv(os.path.join(RES, "real_more_v2.csv"))
    more = more[more.error.isna() & (more.method == "hc")]
    snap = pd.read_csv(os.path.join(RES, "real_detection_v3.csv"))
    snap = snap[(snap.setting == "full") & (snap.method == "hc")]
    rows = []
    for net, g in list(more.groupby("network")) + list(snap.groupby("network")):
        row = dict(network=net)
        for col in ("f1", "f1_ge3"):
            free = g[g.label_free == "max_overlapping_modularity"][col].max()
            best = g[g.label_free.isna()][col].max()
            row.update({f"{col}_free": free, f"{col}_best": best, f"{col}_ratio": free / best})
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "real_label_free_retention.csv"), index=False)


if __name__ == "__main__":
    family_vs_jaccard()
    hetero_by_overlap()
    selftune_vs_oracle()
    real_label_free()
    print("wrote family_vs_jaccard_paired.csv, hetero_by_overlap.csv, selftune_vs_oracle.csv, "
          "real_label_free_retention.csv")
