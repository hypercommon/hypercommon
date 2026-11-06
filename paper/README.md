# Reproducing the paper

Everything below runs from `paper/`. Every number in the manuscript is produced by
the scripts listed here; tables are generated from the result files, never typed.
See `RESULTS_INDEX.md` for where each result lives.

## Software
Python 3.12.10 with numpy 2.4.1, scipy 1.17.1, networkx 3.6.1, pandas 2.3.3,
matplotlib 3.10.8, cdlib 0.4.0, python-igraph 1.0.0, leidenalg 0.11.0, infomap 2.10.1.
The method itself is `../hypercommon/` (repository commit `bcbe7d1` plus the changes
in this folder). LaTeX: MiKTeX 25.12, `revtex4-2`.

## Data
Real networks: SNAP (https://snap.stanford.edu/data/#communities), in `analysis/snap/`:
`com-{dblp,amazon,youtube}.ungraph.txt.gz`, ground truth `com-dblp.all.cmty.txt.gz`,
`com-amazon.all.dedup.cmty.txt.gz`, `com-youtube.all.cmty.txt.gz`, and the top-5000
files `com-{dblp,amazon}.top5000.cmty.txt.gz`. Zachary's karate club from networkx.

The fourteen further networks come from the Netzschleuder repository
(https://networks.skewed.de), CSV format, in `analysis/data_netz/`, each file saved as
`<network>__<subnetwork>.csv.zip` from `https://networks.skewed.de/net/<network>/files/<subnetwork>.csv.zip`:
`ego_social__facebook_{0,107,348,414,686,698,1684,1912,3437,3980}`, `football__football`,
`polbooks__polbooks`, `polblogs__polblogs`, `sp_primary_school__sp_primary_school`.
`analysis/real_more.py` reads them; ground-truth groups of fewer than three nodes are
dropped (and the primary school's teachers are not a class).

LFR graphs: the reference generator of Lancichinetti & Fortunato,
https://github.com/andrealancichinetti/LFRbenchmarks at commit `ec9a860`,
`unweighted_undirected/Sources/benchm.cpp`, compiled unmodified with
`g++ -O3 -std=c++98 -funroll-loops -static` (g++ 16.1, MSYS2) into `bench/bin/`.
The generator reads its seed from `time_seed.dat`; `bench/graphs.py` writes it before
each run, and every generated graph is cached with its arguments in
`bench/lfr_cache/<hash>/args.txt`.

## Seeds
Benchmark graph `(experiment, family, p, seed)` is generated from a string-seeded RNG,
`random.Random("rings|sizes|ks|overlap|p|seed[|rewire]")`, so every graph is rebuilt
exactly. Test seeds are 0–9, calibration seeds 1000–1002. Stochastic methods
(Leiden, Infomap, SLPA, label propagation) are seeded from the graph seed; on the real
networks, Leiden and Infomap use seed 1 and SLPA seeds `random` and `numpy` with 1.

## Timings
Two machines, both Windows 11 (build 26300), Python 3.12.10, NumPy 2.4.1, SciPy 1.17.1.

| Machine | Hardware | Used for |
|---|---|---|
| A | Intel Core Ultra 9 185H (16 cores, 22 threads), 32 GB RAM | everything not listed under B: the benchmarks, the running-time curve (Fig. of Sec. VII I, `bench/runtime.py`) and the single-core re-timing on full Amazon and DBLP (`bench/retime_big.py`) |
| B | workstation, 28 logical processors, about 128 GB RAM | the full SNAP graphs (`bench/big_job.py`, `bench/youtube_job.py`, `bench/slpa_youtube_job.py`) |

Load and threads of each timing:
- `bench/runtime.py`: machine A, run alone; three graphs per size (repetitions as in the
  figure caption).
- `bench/retime_big.py`: machine A, idle, each method in turn with one thread
  (`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `NUMEXPR_NUM_THREADS` = 1).
- `bench/youtube_job.py`: machine B, run alone (log and metadata in
  `bench/results/youtube_full_dl/results/`: 15.5 min for all four HyperCommon variants,
  10.4 min for HyperCommon itself, peak 25.7 GB).
- `bench/big_job.py`: machine B, 16 worker processes at once, so its Leiden and Infomap
  times on YouTube share the machine with up to fifteen other jobs and are indicative only.
- The benchmark runs (`bench/run2.py`, 16 workers) are not timings; their `seconds`
  column is not used.

## Commands
Everything, in order, from scratch: `bash reproduce_all.sh` (about two days on the
computer of the Timings section; the benchmark, detection and timing steps skip work already in their
output files, so an interrupted run can be restarted). The steps it runs:

```
# empirical commonality (Sec. II)
python analysis/real_networks.py            # classes, AUC, degree strata, wedges
python analysis/real_networks_v2.py         # formulas, pair vs triple rule (matched recall)
python analysis/real_family.py              # the commonality family on edges

# synthetic benchmarks: resumable, one CSV per experiment in bench/results/
python bench/run.py <exp>       # exp in rings_uniform, rings_hetero, rings_hetero_overlap, rings_pref,
                                #         lfr, lfr_slices, lfr_clustered, commonality_family
python bench/run.py <ring exp> --methods slpa_x,demon_x,angel_x,cpm_x,scan_x   # competitors' grid extensions
python bench/run.py <lfr exp> --methods slpa100                                # SLPA, 100 iterations
python bench/eval_self_tune.py <exp> 3      # label-free threshold (overlapping modularity)
python bench/analysis.py <exp>              # oracle / frozen protocols, bootstrap CIs, paired differences

# real networks
python bench/real_detection_v3.py <amazon|dblp> <full|induced>   # SNAP, top-5000 ground truth
python bench/real_more_detection_v2.py 16                         # 14 Netzschleuder networks

# running time (run alone on the machine)
python bench/runtime.py

# earlier runs: the comparisons of the earlier tables, kept for reference
python bench/run2.py <exp>      # exp in rings_isolated, rings_uniform, rings_hetero, rings_hetero_overlap,
                                #         rings_pref, lfr, lfr_slices, lfr_clustered, rings_geometry
python bench/analysis2.py <exp> ; python bench/robust_p.py
python bench/run2_family.py ; python bench/run2_family.py --top      # commonality family on the benchmarks
python bench/real_more_detection_v3.py 16                            # 14 Netzschleuder networks
python bench/big_job.py         # SNAP covers, large machine (~100 GB for clique percolation on YouTube)
python bench/youtube_job.py     # HC on full YouTube (machine B of the Timings section)
python bench/slpa_youtube_job.py
python bench/score_big.py ; python bench/score_big.py --onmi-only    # score every SNAP cover
python analysis/make_tables2.py ; python analysis/make_tables2.py family geometry
python analysis/fig_benchmarks2.py

# derived numbers, tables, figures, PDF
python analysis/lfr_triangles.py            # where LFR triangles lie (inside / across communities)
python analysis/derived_numbers.py ; python analysis/export_numbers.py ; python analysis/make_tables.py
python analysis/fig_<name>.py               # every figure script
pdflatex main ; bibtex main ; pdflatex main ; pdflatex main
```

### Main runs (behind the manuscript's tables)
```
# benchmarks: every method, including the native defaults (leiden_g1, leiden_cpm_dens, linkcomm_pd)
python bench/run2.py <exp>      # exp in rings_uniform, rings_hetero, rings_hetero_overlap, rings_pref,
                                #         lfr, lfr_slices, lfr_clustered, rings_geometry, nulls
python bench/run2.py <exp> --methods hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
                                # control families: rings_hubs, rings_split, rings_geometry_swap, rgg
# nested stochastic block model: Linux (WSL) with graph-tool, environment in environment-sbm.yml
#   micromamba create -f environment-sbm.yml ; python bench/run_sbm.py <exp>
#   (same graphs as run2.py; checked edge for edge against the run2 rows)
python bench/analysis2.py <exp> # oracle, label-free, frozen and native protocols, paired bootstrap marks

# theory and implementation checks
python analysis/exact_windows.py            # every recovery window of Sec. V in exact rational arithmetic
python analysis/double_bonds.py             # shared pairs and fused rings (sections/table_doublebonds.tex)
python bench/check_pruning.py               # degree-ratio pruning gives identical covers (Amazon, DBLP)
python bench/retime_big.py                  # single-thread re-timing on Amazon/DBLP; run on an idle machine

# tables and figures
python analysis/make_tables3.py             # the main benchmark tables and the benchmark tables of the Supplemental Material
python analysis/make_tables3_controls.py    # control families and null graphs (Supplemental Material)
python analysis/make_tables3_real.py        # real networks (main table and Supplemental Material)
python analysis/make_tables2.py family geometry
python -c "import sys; sys.path.insert(0, 'analysis'); import make_tables; make_tables.double_bond_table()"
python analysis/fig_<name>.py               # fig_family, fig_real_networks, fig_schematic, fig_omega_t,
                                            # fig_benchmarks2, fig_runtime
```
#### Every table and figure, and the command that makes it
The manuscript is `main.tex`; the Supplemental Material is `supp.tex` (compile
`main` then `supp`, each with `pdflatex`, `bibtex`, `pdflatex`, `pdflatex`; they share labels
through `xr-hyper`). Run every command from `paper/`.

| Item | Command | Output file |
|---|---|---|
| Fig. 1 | `python analysis/real_networks.py ; python analysis/real_networks_v2.py ; python analysis/real_family.py ; python analysis/fig_commonality_combined.py` | `figures/fig_commonality.pdf` |
| Fig. 2 | `python analysis/fig_schematic.py` | `figures/fig_schematic.pdf` |
| Fig. 3 | `python analysis/fig_omega_t.py` | `figures/fig_omega_t.pdf` |
| Fig. 4, Fig. S1 | `python analysis/fig_benchmarks2.py` | `figures/fig_benchmarks.pdf`, `fig_benchmarks_oracle.pdf` |
| Fig. 5 | `python bench/runtime.py ; python analysis/fig_runtime.py` | `figures/fig_runtime.pdf` |
| Figs. S2–S4 | `python analysis/fig_real_networks.py ; python analysis/fig_family.py` | `figures/fig_real_commonality.pdf`, `fig_real_family.pdf`, `fig_real_triples.pdf` |
| Table I (scorecard), Tables S1, S3–S18 | `python analysis/make_tables3.py` | `sections/table_scorecard.tex`, `table_benchmarks.tex`, `supp_benchmarks.tex` |
| Table II (geometry) | `python analysis/make_tables2.py geometry` | `sections/table_geometry.tex` |
| Table III (SNAP), Tables S23–S29 | `python analysis/make_tables3_real.py` | `sections/table_real.tex`, `table_real_induced.tex`, `supp_real.tex` |
| Tables IV, V (fourteen networks) | `python analysis/make_tables2.py more` | `sections/table_real_more.tex` |
| Tables VI, VIII, IX | written by hand from the data sources and `bench/methods2.py` | `sections/appendix.tex` |
| Table VII, Table X | `python analysis/make_tables.py` (`double_bond_table` for Table X, after `python analysis/double_bonds.py`) | `sections/table_data_more.tex`, `table_doublebonds.tex` |
| Table S2 (noise robustness) | `python bench/robust_p.py ; python analysis/make_tables2.py robust` | `sections/table_robust.tex` |
| Tables S19, S20 (controls, null graphs) | `python analysis/make_tables3_controls.py` | `sections/supp_controls.tex` |
| Tables S21, S22 (commonality family) | `python bench/run2_family.py ; python bench/run2_family.py --top ; python analysis/make_tables2.py family` | `sections/table_family.tex` |
| Table S30 (sources of text-only numbers) | written by hand; each row names its script | `sections/supp_diagnostics.tex` |

Which script writes which table: `make_tables3.py` writes `sections/table_benchmarks.tex`
(label-free), `table_scorecard.tex` and `supp_benchmarks.tex` (oracle, frozen, native,
ONMI, all-community recovery, small and large rings, coverage, number of communities,
intervals, selector, grid edges); `make_tables3_controls.py` writes `supp_controls.tex`;
`make_tables3_real.py` writes `table_real.tex` and `supp_real.tex`; `make_tables2.py`
writes `table_family.tex`, `table_geometry.tex`, `table_real_more.tex` and
`table_real_induced.tex`. The numbers behind each table are also written as CSV to
`figures/data/`. The raw result files, one row per graph, method and parameter, are in
the repository compressed: `bench/results2_gz/*.csv.gz` (main runs, behind the manuscript's
tables) and `bench/results_gz/*.csv.gz` (earlier runs and the real networks).
Unpack them into `bench/results2/` and `bench/results/` before running the table
scripts without rerunning the benchmarks; the cached LFR graphs are in
`bench/lfr_cache.tar.xz` (`tar -xJf bench/lfr_cache.tar.xz -C bench`). The SNAP files
are not redistributed; download them as listed under Data.

Scoring on real networks: detected communities are restricted to the nodes that carry
at least one ground-truth group; every nonempty restricted community is kept (primary
rule), and the rule that also drops restricted communities with fewer than three nodes
is recorded alongside (`*_ge3` columns) as a sensitivity check.

## Checks
`python -m pytest ../tests bench/test_hc_fast.py bench/test_hc_scale.py bench/test_hc_scale_family.py`:
the library's tests, plus equality of the fast and scalable implementations (Jaccard and
the whole commonality family) with the reference, and of the triangles-only variant with
Jaccard-filtered 3-clique percolation.
