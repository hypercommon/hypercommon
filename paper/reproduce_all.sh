#!/bin/bash
# Reproduce every result, table, figure and the PDF of the paper, from scratch.
# Run from paper/ with the project's virtual environment active (see README.md).
# Total time on a 22-core desktop: about two days; the benchmark, detection and timing
# steps skip work already in their output files, so an interrupted run can simply be
# restarted. Memory-heavy steps run one at a time.
set -e
cd "$(dirname "$0")"
PY=${PY:-python}
W=${WORKERS:-20}

# --- empirical commonality (Sec. II) --------------------------------------------------
$PY analysis/real_networks.py
$PY analysis/real_networks_v2.py
$PY analysis/real_family.py

# --- synthetic benchmarks (Secs. V-VI) --------------------------------------------------
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered commonality_family; do
  $PY bench/run.py $EXP --workers $W
done
# competitors' grid extensions on the ring families, SLPA with 100 iterations on LFR
for EXP in rings_uniform rings_hetero rings_pref; do
  $PY bench/run.py $EXP --workers $W --methods slpa_x,demon_x,angel_x,cpm_x,scan_x
done
for EXP in lfr lfr_slices lfr_clustered; do
  $PY bench/run.py $EXP --workers 12 --methods slpa100
done
# label-free threshold choice
for EXP in rings_uniform rings_hetero rings_pref lfr; do
  $PY bench/eval_self_tune.py $EXP 3
done
# protocols, bootstrap intervals
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered commonality_family; do
  $PY bench/analysis.py $EXP
done

# --- real networks: detection -------------------------------------------------------------
for NET in amazon dblp; do
  for SETTING in full induced; do
    $PY bench/real_detection_v3.py $NET $SETTING
  done
done
$PY bench/real_more_detection_v2.py 16

# --- running time (alone on the machine) ------------------------------------------------------
$PY bench/runtime.py

# --- earlier runs: the comparisons of the earlier tables ---------------------------
for EXP in rings_isolated rings_hetero rings_hetero_overlap rings_pref rings_uniform lfr lfr_slices lfr_clustered rings_geometry; do
  $PY bench/run2.py $EXP --workers $W
  $PY bench/analysis2.py $EXP
done
$PY bench/robust_p.py
$PY bench/run2_family.py --workers $W
$PY bench/run2_family.py --workers $W --top
$PY bench/real_more_detection_v3.py 16
# SNAP covers: big_job.py needs ~100 GB (clique percolation on YouTube), youtube_job.py ~128 GB;
# copy their outputs to bench/results/big_dl/results_big and bench/results/youtube_full_dl/results
$PY bench/big_job.py
$PY bench/youtube_job.py
$PY bench/slpa_youtube_job.py
$PY bench/score_big.py
$PY bench/score_big.py --onmi-only

# --- derived numbers, tables, figures, PDF ----------------------------------------------------
$PY analysis/lfr_triangles.py
$PY analysis/derived_numbers.py
$PY analysis/export_numbers.py
$PY analysis/make_tables.py
$PY analysis/double_bonds.py
$PY -c "import sys; sys.path.insert(0, 'analysis'); import make_tables; make_tables.double_bond_table()"
# these tables last: they replace the earlier table_benchmarks / table_real_more / table_family
$PY analysis/make_tables2.py
$PY analysis/make_tables2.py family geometry
# --- main runs (behind the manuscript's tables) -----------------------------------------------------------
CTRL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered rings_geometry nulls; do
  $PY bench/run2.py $EXP --workers $W
done
for EXP in rings_hubs rings_split rings_geometry_swap rgg; do
  $PY bench/run2.py $EXP --workers $W --methods $CTRL
done
# nested SBM: needs graph-tool (Linux/WSL, environment-sbm.yml); set SBM_PY to that python
if [ -n "$SBM_PY" ]; then
  for EXP in rings_hetero rings_hetero_overlap rings_uniform rings_pref lfr lfr_slices lfr_clustered nulls rgg; do
    (cd bench && $SBM_PY run_sbm.py $EXP --workers 8)
  done
fi
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered rings_geometry \
           nulls rings_hubs rings_split rings_geometry_swap rgg; do
  (cd bench && $PY analysis2.py $EXP)
done
$PY analysis/exact_windows.py
$PY bench/check_pruning.py
$PY bench/retime_big.py      # alone on the machine
(cd analysis && $PY make_tables3.py && $PY make_tables3_controls.py && $PY make_tables3_real.py)
for F in analysis/fig_*.py; do $PY "$F"; done
# manuscript and Supplemental Material share labels (xr-hyper): build both, then rebuild both
for D in main supp; do pdflatex -interaction=nonstopmode $D && bibtex $D; done
for D in main supp main supp; do pdflatex -interaction=nonstopmode $D; done
echo "done"
