# Where every number lives

Nothing in the paper exists only as a picture. Three layers, from raw to summary:

## 0. Main runs — what the manuscript's comparisons use
The benchmark, real-network and commonality-family comparisons of the paper come from the
main runs below (`bench/run2.py`, results in `bench/results2/`, packed in
`bench/results2_gz/`); the earlier runs (sections 1–3) still supply the empirical
commonality analysis (Sec. III), the running-time curve, the LFR triangle count and the
shared-pair table. `paper/README.md` lists the script behind every table and figure.

| File | Content | Made by |
|---|---|---|
| `bench/results2/<exp>.csv` | one row per graph × method × parameter, with `eq` (overlapping modularity, the label-free criterion) and per-community recovery `rec_mean`, `rec_small`, `rec_large`; exp in `rings_isolated`, `rings_uniform`, `rings_hetero`, `rings_hetero_overlap`, `rings_pref`, `lfr`, `lfr_slices`, `lfr_clustered`, `rings_geometry` | `bench/run2.py` (methods and grids: `bench/methods2.py`) |
| `bench/results2/<exp>_summary2.csv`, `_paired2.csv`, `_frozen_params2.csv` | oracle / frozen / label-free means over the noise range with bootstrap intervals; paired differences from HC | `bench/analysis2.py <exp>` |
| `bench/results2/robust_p.csv` | high-noise mean, Ω at the largest p, p½ per ring family and method | `bench/robust_p.py` |
| `bench/results2/family2.csv`, `family2_top.csv` | HC with each commonality function on the benchmark families, quantile levels 0–0.96 and 0.97–0.995 | `bench/run2_family.py`, `--top` |
| `bench/results/big_dl/results_big/` | covers of every method on SNAP Amazon, DBLP, YouTube (full and induced), produced on a large machine | `bench/big_job.py` (+ `bench/slpa_youtube_job.py` for SLPA on induced YouTube) |
| `bench/results/youtube_full_dl/results/` | HC covers on the full YouTube graph, degree-ratio pruning, machine B of `README.md` (Timings) | `bench/youtube_job.py` (`bench/hc_big.py`) |
| `bench/results/real_big_scores.csv`, `real_big_scores_summary.csv` | every cover above scored: F1, Ω (raw and on distinct truth sets), ONMI on the selected rows, coverage, recovery | `bench/score_big.py` (checked by `bench/test_score_big.py`) |
| `bench/results/real_more_v3.csv` | the fourteen Netzschleuder networks, every method × parameter, with `eq` | `bench/real_more_detection_v3.py` |
| `sections/table_*.tex`, `figures/data/benchmarks_summary2.csv`, `real_more_summary2_{f1,omega}.csv`, `family2_vs_jaccard.csv` | the paper's tables and their numbers | `analysis/make_tables2.py [benchmarks robust real more family geometry]` |
| `figures/fig_benchmarks.pdf`, `fig_benchmarks_oracle.pdf`, `figures/data/benchmarks_curves2.csv` | Ω vs noise, label-free and oracle | `analysis/fig_benchmarks2.py` |

## 1. Raw results (one row per graph × method × parameter)
`bench/results/<experiment>.csv` — experiment, family, p (or μ for LFR), seed, split (test / cal),
method, param, omega, onmi, f1, coverage, n_comms, n_nodes, n_edges, n_wedges, seconds, error.
Experiments: `rings_uniform`, `rings_hetero`, `rings_hetero_overlap` (unequal rings at overlap 0 and
0.10), `rings_pref`, `lfr`, `lfr_slices`, `lfr_clustered`, `commonality_family`. Method names ending
in `_x` are the competitors' grid extensions and `slpa100` is SLPA with 100 iterations; the analysis
merges them into their methods. Ring summaries use p in [0, 0.8].
Every graph is rebuilt exactly from (experiment, family, p, seed) — see `bench/experiments.py`.
LFR graphs are cached as files in `bench/lfr_cache/` (network.dat, community.dat, args.txt).

`bench/results/<exp>_selftune.csv` — the t HyperCommon picks for itself (max overlapping modularity), per graph.

`bench/results/real_detection_v3.csv` — real-network detection used in the paper (SNAP Amazon, DBLP;
truth = top-5000 communities): settings `full` (detect on the whole graph, restrict to the top-5000
nodes; main table) and `induced` (appendix); every method × parameter with best-match F1, ONMI,
coverage of the labelled nodes, n_comms, seconds; every nonempty restricted community kept, `*_ge3`
columns drop restricted communities of fewer than three nodes; label-free rows (self-chosen t for
HyperCommon, γ=1 for Leiden); connected components as a trivial baseline. Data: `analysis/snap/`.
`bench/results/real_more_v2.csv` — 14 further networks (10 Facebook ego networks with circles, football,
political books and blogs, primary school), every method × parameter × seed, F1 / ONMI / Omega /
coverage under both scoring rules; every HyperCommon variant on the same 100 quantile levels; `error`
rows are runs that did not complete.

Superseded, kept for the record, not used in the paper: `real_detection.csv`, `real_detection_v2.csv`
and `real_more.csv` (older scoring: restricted communities of fewer than three nodes always dropped),
`real_family_detection.csv` (an exploratory function sweep on full Amazon/DBLP under that older
scoring), `runtime_tracemalloc_in_timing.csv` (timings taken with the memory tracer active).
`bench/results/runtime.csv` — running times vs n (3 repetitions per size).
`bench/results/slpa_iterations*.csv` — SLPA 21 vs 100 iterations check.

## 2. Real-network analysis
`analysis/out/<Network>.json`, `<Network>_v2.json`, `<Network>_family.json` — all statistics of the
Karate / DBLP / Amazon / YouTube analysis (class counts, quantiles, AUCs, degree strata, pair-vs-triple
precision/recall curves, per-function results).

## 3. Readable tables (what each figure shows)
`figures/data/*.csv`, regenerated by `python analysis/export_numbers.py`:

| File | Content |
|---|---|
| `real_jaccard_by_class.csv` | inner / outer / unlabelled edge counts and Jaccard quantiles; AUCs |
| `real_jaccard_auc_by_degree.csv` | AUC within degree strata |
| `real_family_auc_and_boxes.csv` | every commonality function: AUC (all / stratified) + box quantiles |
| `real_pair_vs_triple_rule.csv` | precision / recall of pair rule vs triple rule per threshold |
| `<exp>_mean_by_param.csv` | mean ± std over test graphs for every method and parameter |
| `<exp>_oracle_by_p.csv` | each method at its best parameter per graph, averaged, per noise level |
| `<exp>_frozen_params.csv` / `<exp>_frozen_by_p.csv` | parameter chosen on calibration graphs, and its test scores |
| `rings_uniform_selftuned_by_p.csv` | HyperCommon with its self-chosen t |
| `selftune_vs_oracle.csv` | label-free t vs oracle t (mean Ω, ratio) and the best other method under the frozen protocol, every family |
| `real_more_summary.csv` | 14 Netzschleuder networks: HC, best HC_f, best other, under both scoring rules |
| `real_label_free_retention.csv` | 16 real networks: HC F1 at the label-free threshold over F1 at the best threshold, both rules |
| `runtime_summary.csv` / `runtime_exponents.csv` | running time per size (mean of 3) and fitted exponents |
| `lfr_triangles.csv` | per LFR test graph: triangles, share with all three nodes in one planted community, transitivity |

## Code that produced them
`bench/` (harness, methods, metrics, graphs), `analysis/` (real networks, figures, export).
Tests: `bench/test_hc_fast.py` (fast implementation = library; triangles-only = filtered 3-CPM).
