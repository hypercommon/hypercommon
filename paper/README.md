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

Wedge counts (Table S-wedges, `analysis/wedge_scaling.py`): 29 further SNAP edge lists,
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
Every result behind the manuscript and the Supplemental Material was computed on one
machine, a Mac mini (Apple M4 Pro, 14 cores, 64 GB RAM, macOS 27.2), in three micromamba
environments under `~/hc/mamba`: `hc` (Python 3.12.16, numpy 2.4.1, scipy 1.17.1,
networkx 3.6.1, python-igraph 1.0.0, leidenalg 0.11.0, infomap 2.10.1, cdlib 0.4.0) for
everything except the two below; `gt` (graph-tool 3.9) for the nested stochastic block
model; `cpm` (numba 0.68.1) for `bench/cpm_bits.py`. OSLOM is the OSLOM2 binary built by
`bench/build_oslom.sh`. The manuscript reports no running times; the `seconds` columns of
the result files are kept for reference only. Against an earlier run of the same
scripts under Windows and Linux, every deterministic method gives identical scores;
Infomap's random numbers depend on the platform (mean Omega 0.831 against 0.837), and
the stochastic block model is stochastic.

`bench/queue_mac.sh` runs the whole sequence on that machine (resumable, one done marker
per step), `bench/native_mac.sh` the native-default runs beside it, and
`bench/watchdog_mac.sh` kills a run above 30 GB of memory or when swap exceeds 12 GB
(an OSLOM process alone, so that its Python caller records the failure).

## Commands
Everything, in order, from scratch: `bash reproduce_all.sh` (several days on the
machine above; the benchmark and detection steps skip work already in their
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
| Figs. S2–S4 | `python analysis/fig_real_networks.py ; python analysis/fig_family.py` | `figures/fig_real_commonality.pdf`, `fig_real_family.pdf`, `fig_real_triples.pdf` |
| Table I (scorecard), Tables S1, S3–S18 | `python analysis/make_tables3.py` | `sections/table_scorecard.tex`, `table_benchmarks.tex`, `supp_benchmarks.tex` |
| Table II (geometry) | `python analysis/make_tables2.py geometry` | `sections/table_geometry.tex` |
| Table III (SNAP), Tables S23, S25–S30 | `python analysis/make_tables3_real.py` | `sections/table_real.tex`, `table_real_induced.tex`, `supp_real.tex` |
| Tables IV, V (fourteen networks) | `python analysis/make_tables2.py more` | `sections/table_real_more.tex` |
| Tables VI, VIII, IX | written by hand from the data sources and `bench/methods2.py` | `sections/appendix.tex` |
| Table VII, Table X | `python analysis/make_tables.py` (`double_bond_table` for Table X, after `python analysis/double_bonds.py`) | `sections/table_data_more.tex`, `table_doublebonds.tex` |
| Table S2 (noise robustness) | `python bench/robust_p.py ; python analysis/make_tables2.py robust` | `sections/table_robust.tex` |
| Tables S19, S20 (controls, null graphs) | `python analysis/make_tables3_controls.py` | `sections/supp_controls.tex` |
| Tables S21, S22 (commonality family) | `python bench/run2_family.py ; python bench/run2_family.py --top ; python analysis/make_tables2.py family` | `sections/table_family.tex` |
| Table S24 (wedges of 32 SNAP graphs) | `python analysis/wedge_scaling.py` | `sections/table_wedges.tex` |
| Table S31 (sources of text-only numbers) | written by hand; each row names its script | `sections/supp_diagnostics.tex` |

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
