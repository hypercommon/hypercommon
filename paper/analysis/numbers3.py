"""
Numbers quoted in the text of the third version that are not read directly off one
table: counts of settings by paired mark for every metric and protocol, the label-free
to oracle ratio, the ablation counts. Writes figures/data/numbers3.csv and prints it.

    python numbers3.py
"""
import os

import pandas as pd

from make_tables3 import DATA, RING_ROWS, ROWS, load, verdict

S, P = load()
out = []


def put(key, value, note=""):
    out.append(dict(key=key, value=value, note=note))


# marks against the best other method, per metric and protocol
for metric in ("omega", "onmi", "rec_mean"):
    for prot in ("oracle", "label_free", "frozen", "native"):
        for name, rows in (("rings", RING_ROWS), ("lfr", ROWS[len(RING_ROWS):])):
            V = [verdict(S, P, e, f, prot, metric) for e, f, *_ in rows]
            V = [v for v in V if v]
            for m in "+=-":
                put(f"{metric}.{prot}.{name}.{m}", sum(v["mark"] == m for v in V), f"of {len(V)}")

# label-free threshold against the best threshold, HC
r = []
for e, f, *_ in ROWS:
    a, b = S[(S.experiment == e) & (S.family == f) & (S.metric == "omega") & (S.method == "hc")
             & S.protocol.isin(["oracle", "label_free"])].set_index("protocol")["mean"][["label_free", "oracle"]]
    r.append(dict(family=f, ratio=a / b))
r = pd.DataFrame(r)
put("ratio.n", len(r))
put("ratio.ge95", int((r.ratio >= 0.95).sum()))
put("ratio.ge90", int((r.ratio >= 0.90).sum()))
put("ratio.min", round(r.ratio.min(), 3), r.loc[r.ratio.idxmin(), "family"])

# ablations: HC against each ablation, paired, every benchmark setting
for abl in ("hc_2test", "hc_tri", "hc_econn"):
    for prot in ("oracle", "label_free", "frozen"):
        for name, rows in (("rings", RING_ROWS), ("lfr", ROWS[len(RING_ROWS):]), ("all", ROWS)):
            d = []
            for e, f, *_ in rows:
                p = P[(P.experiment == e) & (P.family == f) & (P.protocol == prot) & (P.metric == "omega")
                      & (P.method == abl)]
                if len(p):
                    d.append(p.iloc[0])
            d = pd.DataFrame(d)
            # diff = ablation - HC
            put(f"abl.{abl}.{prot}.{name}.hc_better", int((d.hi < 0).sum()), f"of {len(d)}")
            put(f"abl.{abl}.{prot}.{name}.abl_better", int((d.lo > 0).sum()), f"of {len(d)}")
            put(f"abl.{abl}.{prot}.{name}.abl_higher_mean", int((d["diff"] > 0).sum()), f"of {len(d)}")
            put(f"abl.{abl}.{prot}.{name}.loss_range", f"{-d['diff'].max():.2f}..{-d['diff'].min():.2f}")
# link communities against HC
for prot in ("oracle", "label_free"):
    for name, rows in (("rings", RING_ROWS), ("lfr", ROWS[len(RING_ROWS):])):
        d = []
        for e, f, *_ in rows:
            p = P[(P.experiment == e) & (P.family == f) & (P.protocol == prot) & (P.metric == "omega")
                  & (P.method == "linkcomm")]
            if len(p):
                d.append(dict(family=f, diff=p.iloc[0]["diff"], lo=p.iloc[0].lo, hi=p.iloc[0].hi))
        d = pd.DataFrame(d)
        put(f"lk.{prot}.{name}.hc_better", int((d.hi < 0).sum()), f"of {len(d)}")
        put(f"lk.{prot}.{name}.lk_better", ";".join(d[d.lo > 0].family), "")

D = pd.DataFrame(out)
D.to_csv(os.path.join(DATA, "numbers3.csv"), index=False)
pd.set_option("display.max_rows", 500, "display.width", 200)
print(D.to_string(index=False))
