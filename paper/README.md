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

Lattice-like real networks (Sec. VII G, Table II, Table S26): eleven kinds, each built
from public data by `analysis/<kind>_build.py` into `analysis/data_<kind>/`, one folder
per graph (`edges.tsv`, `groups.txt`) and its position-only control (`<name>_null/`):
`hic` (Rao et al. 2014, GEO GSE63525; `hic_null.py` writes the controls), `protein`
(CATH v4.4 domains on PDB chains; `protein_null.py`), `sketch` (Quick, Draw!), `image`
(BSDS500), `tracks` (MOT17 labels), `geo` (US ACS commuting flows and UK ONS ODWP01EW),
`transit` (twelve agency GTFS feeds), `streets` (OpenStreetMap via Overpass), `tssb`
(Time Series Segmentation Benchmark), `spatial` (DLPFC Visium slices, Maynard et al.
2021) and `text` (English Wikipedia featured articles). The built graphs are not in the
repository, because the source licences differ (BSDS500 is research-only, MOT17 is
CC BY-NC-SA, OpenStreetMap and some feeds are share-alike, several transit feeds allow
no redistribution). Each folder's `README.md` is kept: the exact URLs, download dates,
checksums, licence and attribution of its sources, and the build command. Rebuild with
`python -I analysis/<kind>_build.py ...` as given there.

Wedge counts (Table S21, `analysis/wedge_scaling.py`): 29 further SNAP edge lists,
no ground truth, in `analysis/snap/more/`, each saved under its SNAP file name from
`https://snap.stanford.edu/data/<file>`: `ca-{GrQc,HepTh,HepPh,AstroPh,CondMat}`,
`email-{Enron,EuAll}`, `wiki-{Talk,Vote}`, `facebook_combined`, `twitter_combined`,
`soc-Epinions1`, `soc-Slashdot0902`, `loc-{brightkite,gowalla}_edges`, `com-lj.ungraph`
(under `bigdata/communities/`), `as-skitter`, `as-caida20071105`, `oregon1_010526`,
`p2p-Gnutella31`, `web-{Google,Stanford,NotreDame,BerkStan}`, `roadNet-{CA,PA}`,
`cit-{HepPh,Patents}`, `amazon0601`, all `.txt.gz`; the counts are in
`figures/data/wedge_scaling.csv`.

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

## Machine
Every result behind the manuscript and the Supplemental Material, except the later parts
named at the end of this section, was computed on one machine, a Mac mini (Apple M4 Pro, 14 cores, 64 GB RAM, macOS 27.2), in three micromamba
environments under `~/hc/mamba`: `hc` (Python 3.12.16, numpy 2.4.1, scipy 1.17.1,
networkx 3.6.1, python-igraph 1.0.0, leidenalg 0.11.0, infomap 2.10.1, cdlib 0.4.0) for
everything except the two below; `gt` (graph-tool 3.9) for the nested stochastic block
model; `cpm` (numba 0.68.1) for `bench/cpm_bits.py`. OSLOM is the OSLOM2 binary built by
`bench/build_oslom.sh`. The `seconds` columns of the result files are kept for reference
only (runs in parallel). The manuscript reports no running times; controlled timings
(Running times, under Commands) are kept in the repository, not in the paper. Against an earlier run of the same
scripts under Windows and Linux, every deterministic method gives identical scores;
Infomap's random numbers depend on the platform (mean Omega 0.831 against 0.837), and
the stochastic block model is stochastic.

`bench/queue_mac.sh` runs the whole sequence on that machine (resumable, one done marker
per step), `bench/native_mac.sh` the native-default runs beside it, and
`bench/watchdog_mac.sh` kills a run above 30 GB of memory or when swap exceeds 12 GB
(an OSLOM process alone, so that its Python caller records the failure).

Two later parts (and the two methods added last, below) were computed on a second machine, a Windows 11 laptop (Intel Core Ultra 9
185H, 22 threads, 32 GB RAM): the robustness sweep (`bench/queue_dense.sh`, Python 3.12.10
with the versions under Software) and the lattice-like real networks
(`bench/real_lattice_detection.py`); OSLOM and the SBM for both ran under WSL (Ubuntu
24.04, Python 3.12.16, graph-tool 3.9, OSLOM2 built by `bench/build_oslom.sh`). The
robustness curves join the main runs at p = 0, 0.05, ..., 0.8 (Mac) with the dense points
(laptop); the deterministic methods agree between the two machines (see above).

BigClam and ego-splitting, the fifteenth and sixteenth methods compared, were added last,
on the laptop, for every part of the comparison (Python 3.12.10 with numba 0.68.1 and
llvmlite 0.50.0 beside the versions under Software): `bench/bigclam.py` is a port of SNAP's
`agmfast` (BigClam), `bench/egosplit.py` of the ego-splitting reference code.

## Commands
Everything, in order, from scratch: `bash reproduce_all.sh` (several days on the
machine above; the benchmark and detection steps skip work already in their
output files, so an interrupted run can be restarted). The steps it runs:

```
# empirical commonality (Sec. III)
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

# earlier runs: the comparisons of the earlier tables, kept for reference
python bench/run2.py <exp>      # exp in rings_isolated, rings_uniform, rings_hetero, rings_hetero_overlap,
                                #         rings_pref, lfr, lfr_slices, lfr_clustered, rings_geometry
python bench/analysis2.py <exp> ; python bench/robust_p.py
python bench/run2_family.py ; python bench/run2_family.py --top      # commonality family on the benchmarks
python bench/real_more_detection_v3.py 16                            # 14 Netzschleuder networks
python bench/big_job.py         # SNAP covers; clique percolation through bench/cpm_fast.py (bounded memory)
python bench/youtube_job.py     # HC on full YouTube
python bench/slpa_youtube_job.py
python bench/score_big.py ; python bench/score_big.py --onmi-only    # score every SNAP cover
python analysis/make_tables2.py ; python analysis/make_tables2.py family geometry
python analysis/fig_benchmarks2.py

# derived numbers, tables, figures, PDF
python analysis/lfr_triangles.py            # where LFR triangles lie (inside / across communities)
python analysis/wedge_scaling.py            # W against m ln m on 32 SNAP graphs (--table-only: from the CSV)
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
# OSLOM: Linux or macOS; bash bench/build_oslom.sh downloads OSLOM2 from oslom.org and compiles it in ~/tools
#   bash bench/queue_oslom.sh       # python bench/run_oslom.py <exp> over the benchmark experiments
#                                   # (rows in results2/<exp>_oslom.csv)
#   python bench/real_more_detection_v3.py 10   # adds OSLOM (five seeds) to the 14 Netzschleuder networks
#   python bench/run_oslom_big.py --workers 6 --hours 48   # SNAP full and induced graphs, seed 1
#   python bench/score_big.py ; python bench/score_big.py --onmi-only
# clique percolation on Facebook ego 107 (more than 10^7 maximal cliques): numba, environment cpm
#   NUMBA_NUM_THREADS=12 python bench/cpm_bits.py facebook_107 3 4 5 6 7 --save covers.pkl
#   python bench/cpm_rows_from_bits.py facebook_107 covers.pkl   # its rows in real_more_v2/v3
# BigClam and ego-splitting (added last; resumable)
#   bash bench/queue_more.sh            # every benchmark experiment, then the robustness sweep at step 0.05
#   bash bench/queue_more_lattice.sh    # the lattice-like real networks
#   python bench/real_more_detection_v3.py 16    # the Netzschleuder networks (picks up the missing methods)
#   python bench/run_more_big.py --workers 4     # SNAP Amazon, DBLP, YouTube; then score_big.py as above
python bench/analysis2.py <exp> # oracle, label-free, frozen and native protocols, paired bootstrap marks

# theory and implementation checks
python analysis/exact_windows.py            # every recovery window of Sec. V in exact rational arithmetic
python analysis/double_bonds.py             # shared pairs and fused rings (sections/table_doublebonds.tex)
python bench/check_pruning.py               # degree-ratio pruning gives identical covers (Amazon, DBLP)

# tables and figures
python analysis/make_tables3.py             # the main benchmark tables and the benchmark tables of the Supplemental Material
python analysis/make_tables3_controls.py    # control families and null graphs (Supplemental Material)
python analysis/make_tables3_real.py        # real networks (main table and Supplemental Material)
python analysis/make_tables2.py family geometry
python -c "import sys; sys.path.insert(0, 'analysis'); import make_tables; make_tables.double_bond_table()"
python analysis/fig_<name>.py               # fig_family, fig_real_networks, fig_schematic, fig_omega_t,
                                            # fig_benchmarks2
```
### Robustness sweep (Sec. VII C)
```
# p = 0, 0.01, ..., 1 on the test graphs of the 28 ring settings; resumable
bash bench/queue_dense.sh          # every method but OSLOM and the SBM (DENSE_P=1 python bench/run2.py <exp>),
                                   # then python bench/robust_structure.py (p_weak, p_clust)
bash bench/queue_dense_wsl.sh      # OSLOM and the SBM, Linux, p in steps of 0.05 (DENSE_STEP)
python analysis/robust_dense.py    # curves, p_0.6, area, markers; Fig. 4, Figs. S6, S7
python analysis/robust_table.py    # Table S19, paired bootstrap over the test graphs
```
The rows these scripts read (only the columns they use) are in
`bench/results2_gz/<exp>_dense*.csv.gz` and `robust_structure.csv.gz`; unpack them into
`bench/results2/` to rerun the analysis without the sweep.

### Lattice-like real networks (Sec. VII G)
```
python -I analysis/<kind>_build.py ...                 # see analysis/data_<kind>/README.md
python analysis/hic_null.py ; python analysis/protein_null.py   # controls (other kinds: in the builder)
python bench/real_lattice_detection.py <kind> 16 all   # every method but OSLOM; results/real_lattice_<kind>.csv
RL_SUFFIX=_wsl python bench/real_lattice_detection.py <kind> 4 oslom,sbm   # Linux: OSLOM, SBM
python analysis/real_lattice_summary.py <kind>         # per-kind summary, figures/data/real_lattice_summary_<kind>.csv
python analysis/real_lattice_ranks.py                  # Table II, Table S26
python bench/edge_prune.py analysis/snap/com-dblp.ungraph.txt.gz analysis/snap/com-amazon.ungraph.txt.gz \
    analysis/data_spatial/dlpfc_151507_k16/edges.tsv   # edge-pruning numbers of Sec. VIII, figures/data/edge_prune.csv
```

The detection rows (one per graph, method and parameter; scores only, no graph data) are
in `bench/results_gz/real_lattice_<kind>{,_wsl}.csv.gz`; unpack them into `bench/results/`
to rebuild the tables without the graphs.

#### Every table and figure, and the command that makes it
The manuscript is `main.tex`; the Supplemental Material is `supp.tex` (compile
`main` then `supp`, each with `pdflatex`, `bibtex`, `pdflatex`, `pdflatex`; they share labels
through `xr-hyper`). Run every command from `paper/`.

| Item | Command | Output file |
|---|---|---|
| Fig. 1 | `python analysis/real_networks.py ; python analysis/real_networks_v2.py ; python analysis/real_family.py ; python analysis/fig_commonality_combined.py` | `figures/fig_commonality.pdf` |
| Fig. 2 | `python analysis/fig_schematic.py` | `figures/fig_schematic.pdf` |
| Fig. 3 | `python analysis/fig_omega_t.py` | `figures/fig_omega_t.pdf` |
| Fig. 4, Figs. S6, S7 | robustness sweep (below) ; `python analysis/robust_dense.py` | `figures/fig_robust.pdf`, `fig_robust_all.pdf`, `fig_robust_oracle.pdf` |
| Table I (geometry) | `python analysis/make_tables2.py geometry` | `sections/table_geometry.tex` |
| Table II (lattice-like real networks), Table S26 | lattice-like real networks (below) ; `python analysis/real_lattice_ranks.py` | `sections/table_real_lattice.tex`, `table_real_lattice_sets.tex` |
| Figs. S1–S3 | `python analysis/fig_real_networks.py ; python analysis/fig_family.py` | `figures/fig_real_commonality.pdf`, `fig_real_family.pdf`, `fig_real_triples.pdf` |
| Figs. S4, S5 | `python analysis/fig_benchmarks2.py` | `figures/fig_benchmarks.pdf`, `fig_benchmarks_oracle.pdf` |
| Tables S1, S3, S4 (data, ring families, grids) | written by hand from the data sources and `bench/methods2.py` | `sections/appendix.tex` |
| Table S2 (fourteen networks) | `python analysis/make_tables.py` | `sections/table_data_more.tex` |
| Table S5 (shared pairs) | `python analysis/double_bonds.py ; python -c "import sys; sys.path.insert(0, 'analysis'); import make_tables; make_tables.double_bond_table()"` | `sections/table_doublebonds.tex` |
| Table S6 (label-free), Table S7 (scorecard), Tables S8–S14 | `python analysis/make_tables3.py` | `sections/table_benchmarks.tex`, `table_scorecard.tex`, `supp_benchmarks.tex` |
| Tables S15, S16 (controls, null graphs) | `python analysis/make_tables3_controls.py` | `sections/supp_controls.tex` |
| Tables S17, S18 (commonality family) | `python bench/run2_family.py ; python bench/run2_family.py --top ; python analysis/make_tables2.py family` | `sections/table_family.tex` |
| Table S19 (robustness) | after `robust_dense.py`: `python analysis/robust_table.py` | `sections/table_robust_dense.tex` |
| Table S20 (SNAP), Tables S24, S25 | `python analysis/make_tables3_real.py` | `sections/table_real.tex`, `supp_real.tex` |
| Table S21 (wedges of 32 SNAP graphs) | `python analysis/wedge_scaling.py` | `sections/table_wedges.tex` |
| Tables S22, S23 (fourteen networks) | `python analysis/make_tables2.py more` | `sections/table_real_more.tex`, `table_real_more_omega.tex` |
| Table S27 (per-branch thresholds) | runs below (Per-branch thresholds); `python bench/make_table_tree.py` | `sections/table_tree.tex` |

Further tables kept in the repository but not in the paper: `sections/table_times.tex`
(running times, `bench/make_table_times.py`), `sections/extra_benchmarks.tex`
(frozen and native protocols, `make_tables3.py`), `sections/table_real_induced.tex`
(SNAP graphs on the subgraphs induced by the labelled nodes, `make_tables2.py real`) and
`sections/supp_diagnostics.tex` (where numbers quoted in the text without a table come from).

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

Running times (not in the paper; laptop, idle, one thread): `python bench/time_methods.py`
(Windows: every method but OSLOM and the SBM) and, under WSL with the `gt` environment,
`python bench/time_methods.py oslom sbm`, rows to `bench/results/time_methods.csv`;
`RETIME_OUT=retime_big_laptop.csv python bench/retime_big.py amazon dblp` (full SNAP graphs);
then `python bench/make_table_times.py` writes `sections/table_times.tex`.

Per-branch thresholds (Discussion, "a first try"; Table S27): `bench/hc_tree.py`.
`bash bench/run_tree_all.sh` (rings and LFR, rows to `bench/results_tree/`) and
`bash bench/run_tree_real.sh` (the real networks), then `python bench/compare_tree.py` and
`python bench/compare_tree_real.py` write `bench/results_tree/compare_bench.csv` (the
label-free LFR gain of 0.11 is the mean of `tree_minus_hc` over the eight LFR families)
and `compare_real.csv`; `python bench/make_table_tree.py` writes `sections/table_tree.tex`
(Table S27). The raw rows are in `bench/results_tree_gz/*.csv.gz`.

## Checks
`python -m pytest ../tests bench/test_hc_fast.py bench/test_hc_scale.py bench/test_hc_scale_family.py`:
the library's tests, plus equality of the fast and scalable implementations (Jaccard and
the whole commonality family) with the reference, and of the triangles-only variant with
Jaccard-filtered 3-clique percolation.
