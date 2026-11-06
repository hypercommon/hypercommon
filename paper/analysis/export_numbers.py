"""Write every number behind the paper's figures and tables to plain CSV files in
figures/data/, so results can be read without opening any plot. Re-run after any
new result: python analysis/export_numbers.py"""
import json, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
BENCH = os.path.join(HERE, "..", "bench", "results")
DATA = os.path.join(HERE, "..", "figures", "data")
os.makedirs(DATA, exist_ok=True)
NETS = ["Karate", "Dblp", "Amazon", "Youtube"]

# real networks: classes, medians, AUCs, degree strata
rows, strata = [], []
for n in NETS:
    d = json.load(open(os.path.join(OUT, f"{n}.json")))
    e = d["edges"]
    for cls in ("inner", "outer", "unlabelled"):
        q = e["jaccard"].get(cls)
        rows.append(dict(network=n, edge_class=cls, n_edges=e["class_counts"][cls],
                         **({k: q[k] for k in ("p05", "q1", "median", "q3", "p95", "mean")} if q else {})))
    for s in e["degree_stratified"]:
        strata.append(dict(network=n, **s))
    rows.append(dict(network=n, edge_class="AUC_inner_vs_outer", median=e["auc_inner_vs_outer"]))
    rows.append(dict(network=n, edge_class="AUC_inner_vs_outer_plus_unlabelled", median=e["auc_inner_vs_pooled"]))
pd.DataFrame(rows).to_csv(os.path.join(DATA, "real_jaccard_by_class.csv"), index=False)
pd.DataFrame(strata).to_csv(os.path.join(DATA, "real_jaccard_auc_by_degree.csv"), index=False)

# real networks: commonality family
fam = []
for n in NETS:
    d = json.load(open(os.path.join(OUT, f"{n}_family.json")))
    for f, x in d["f"].items():
        fam.append(dict(network=n, function=f, auc_all_edges=x["auc"], auc_degree_stratified=x["auc_strat_mean"],
                        **{f"inner_{k}": v for k, v in zip(["p05", "q1", "median", "q3", "p95"], x["box"]["inner"])},
                        **{f"outer_{k}": v for k, v in zip(["p05", "q1", "median", "q3", "p95"], x["box"]["outer"])}))
pd.DataFrame(fam).to_csv(os.path.join(DATA, "real_family_auc_and_boxes.csv"), index=False)

# real networks: pair rule vs triple rule (unique-triple weighted)
pr = []
for n in NETS:
    w = json.load(open(os.path.join(OUT, f"{n}_v2.json")))["wedges"]
    for rule in ("pair", "triple"):
        for t, x in zip(w["t"], w[rule]):
            pr.append(dict(network=n, rule=rule, t=t, **x))
    pr.append(dict(network=n, rule="triangles_only", **w["triangle_only"]))
    pr.append(dict(network=n, rule="base_rate", precision=w["base_rate"]))
pd.DataFrame(pr).to_csv(os.path.join(DATA, "real_pair_vs_triple_rule.csv"), index=False)

# benchmarks: Omega vs t curves (HyperCommon) and per-experiment summaries
for exp in ("rings_uniform", "rings_hetero", "lfr", "lfr_slices", "lfr_clustered", "rings_pref", "commonality_family"):
    p = os.path.join(BENCH, f"{exp}.csv")
    if not os.path.exists(p):
        continue
    df = pd.read_csv(p)
    df = df[df.error.isna()].copy()
    if exp.startswith("rings"):
        # the paper's noise range p in [0, 0.8], applied before the frozen calibration as in
        # bench/analysis.py (rings_pref also sampled p up to 1.0; those rows stay in the raw file)
        df = df[df.p <= 0.8 + 1e-9]
    df["method"] = df["method"].str.replace(r"_x$", "", regex=True)   # grid extensions
    s100 = df["method"] == "slpa100"                                    # SLPA, 100 iterations
    df.loc[s100, "param"] = df.loc[s100, "param"].astype(str) + "|t100"
    df.loc[s100, "method"] = "slpa"
    mets = [m for m in ("omega", "onmi", "f1", "coverage") if m in df.columns]
    # every (family, method, param, p) averaged over test graphs
    agg = df[df.split == "test"].groupby(["family", "method", "param", "p"])[mets + ["n_comms", "seconds"]].agg(["mean", "std", "count"])
    agg.columns = ["_".join(c) for c in agg.columns]
    agg.reset_index().to_csv(os.path.join(DATA, f"{exp}_mean_by_param.csv"), index=False)
    # oracle (best param per test graph) and frozen (calibration-chosen param)
    test = df[df.split == "test"]
    orc = test.loc[test.groupby(["family", "method", "p", "seed"]).omega.idxmax()]
    orc.groupby(["family", "method", "p"])[mets].mean().reset_index().to_csv(
        os.path.join(DATA, f"{exp}_oracle_by_p.csv"), index=False)
    cal = df[df.split == "cal"].groupby(["family", "method", "param"]).omega.mean().reset_index()
    best = cal.loc[cal.groupby(["family", "method"]).omega.idxmax()][["family", "method", "param"]]
    best.to_csv(os.path.join(DATA, f"{exp}_frozen_params.csv"), index=False)
    fro = test.merge(best, on=["family", "method", "param"])
    fro.groupby(["family", "method", "p"])[mets].mean().reset_index().to_csv(
        os.path.join(DATA, f"{exp}_frozen_by_p.csv"), index=False)

# self-tuned thresholds
st = os.path.join(BENCH, "rings_uniform_selftune.csv")
if os.path.exists(st):
    s = pd.read_csv(st); s["t"] = s.t.round(2); s["p"] = s.p.round(2)
    df = pd.read_csv(os.path.join(BENCH, "rings_uniform.csv")); df = df[(df.split == "test") & df.error.isna()]
    df["t"] = pd.to_numeric(df.param, errors="coerce").round(2); df["p"] = df.p.round(2)
    m = df.merge(s.rename(columns={"t": "t_selected"}), left_on=["family", "method", "p", "seed", "t"],
                 right_on=["family", "method", "p", "seed", "t_selected"])
    m.groupby(["family", "method", "criterion", "p"]).omega.mean().reset_index().to_csv(
        os.path.join(DATA, "rings_uniform_selftuned_by_p.csv"), index=False)
print("wrote", sorted(os.listdir(DATA)))
