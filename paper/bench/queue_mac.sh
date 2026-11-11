#!/bin/bash
# Every run of the paper on one machine (the Mac mini, M4 Pro, 64 GB), largest first.
# Timing steps run alone; the rest run in order. Each step leaves runs/done/<step>, so a
# restarted queue skips finished steps. Usage, from paper/bench:
#   PY=~/hc/mamba/envs/hc/bin/python GT_PY=~/hc/mamba/envs/gt/bin/python bash queue_mac.sh
cd "$(dirname "$0")"
PY=${PY:-python}
W=${WORKERS:-13}
RUNS=${RUNS:-$HOME/hc/runs}
mkdir -p "$RUNS/done" "$RUNS/logs" results results2
export BIG_MEM_LIMIT_GB=${BIG_MEM_LIMIT_GB:-48}
ONE="OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1"

step() {   # step <name> <command...>
  local name=$1; shift
  [ -e "$RUNS/done/$name" ] && { echo "skip $name"; return; }
  echo "$(date '+%F %T') start $name" >> "$RUNS/queue.log"
  if (eval "$@") >> "$RUNS/logs/$name.log" 2>&1; then
    touch "$RUNS/done/$name"; echo "$(date '+%F %T') done  $name" >> "$RUNS/queue.log"
  else
    echo "$(date '+%F %T') FAIL  $name (exit $?)" >> "$RUNS/queue.log"
  fi
}

# --- timings, alone on the machine ------------------------------------------------------
step youtube_hc  "env $ONE $PY -u youtube_job.py --out results/youtube_full_dl/results"
step retime_big  "env $ONE $PY -u retime_big.py amazon dblp youtube"
step runtime     "env $ONE $PY -u runtime.py"

# --- SNAP graphs: covers of every method, then scoring ------------------------------------
step big_job     "$PY -u big_job.py --workers 4 --out results/big_dl/results_big"
step slpa_youtube "$PY -u slpa_youtube_job.py"
step score_big   "$PY -u score_big.py"
step score_onmi  "$PY -u score_big.py --onmi-only"
step check_pruning "$PY -u check_pruning.py"

# --- real networks ---------------------------------------------------------------------------
step real_networks    "cd ../analysis && $PY -u real_networks.py"
step real_networks_v2 "cd ../analysis && $PY -u real_networks_v2.py"
step real_family      "cd ../analysis && $PY -u real_family.py"
for NET in amazon dblp; do for S in full induced; do
  step real_det_${NET}_$S "$PY -u real_detection_v3.py $NET $S"
done; done
step real_more_v2 "$PY -u real_more_detection_v2.py $W"
step real_more_v3 "$PY -u real_more_detection_v3.py $W"
# clique percolation on facebook_107 (2.2 million maximal cliques): bitset version in the
# numba env (CPM_PY), then scored into real_more_v2.csv and real_more_v3.csv
if [ -n "$CPM_PY" ]; then
  step cpm_fb107 "NUMBA_NUM_THREADS=12 $CPM_PY -u cpm_bits.py facebook_107 3 4 5 6 7 --save $RUNS/cpm107/covers.pkl && $PY -u cpm_rows_from_bits.py facebook_107 $RUNS/cpm107/covers.pkl"
fi

# --- synthetic benchmarks: main runs ---------------------------------------------------------------
for EXP in rings_hetero_overlap rings_uniform rings_hetero rings_pref lfr lfr_slices lfr_clustered rings_geometry \
           rings_isolated nulls; do
  step run2_$EXP "$PY -u run2.py $EXP --workers $W"
done
CTRL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
for EXP in rings_hubs rings_split rings_geometry_swap rgg; do
  step run2_$EXP "$PY -u run2.py $EXP --workers $W --methods $CTRL"
done
# native defaults (Leiden at gamma = 1, Leiden-CPM at the edge density, link communities at
# the cut of largest partition density) on the main experiments; DEMON on the swapped geometry
NATIVE=leiden_g1,leiden_cpm_dens,linkcomm_pd
for EXP in rings_hetero_overlap rings_uniform rings_hetero rings_pref lfr lfr_slices lfr_clustered rings_geometry nulls; do
  step run2native_$EXP "$PY -u run2.py $EXP --workers $W --methods $NATIVE"
done
step run2_rings_geometry_swap_demon "$PY -u run2.py rings_geometry_swap --workers $W --methods demon"
step robust_p    "$PY -u robust_p.py"
step run2_family "$PY -u run2_family.py --workers $W"
step run2_family_top "$PY -u run2_family.py --workers $W --top"

# --- stochastic block model and OSLOM (graph-tool env; OSLOM2 binary) ----------------------------
if [ -n "$GT_PY" ]; then
  for EXP in rings_hetero_overlap rings_hetero rings_uniform rings_pref lfr lfr_slices lfr_clustered nulls rgg; do
    step sbm_$EXP "$GT_PY -u run_sbm.py $EXP --workers $W"
  done
fi
if [ -x "${OSLOM_BIN:-$HOME/tools/OSLOM2/oslom_undir}" ]; then
  for EXP in lfr rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr_slices lfr_clustered nulls rgg; do
    step oslom_$EXP "$PY -u run_oslom.py $EXP --workers $W"
  done
fi
# OSLOM on the real networks: the fourteen smaller ones (added to real_more_v3.csv), then
# Amazon, DBLP and YouTube (big_job.py covers, at most 48 h per run), scored with the rest
if [ -x "${OSLOM_BIN:-$HOME/tools/OSLOM2/oslom_undir}" ]; then
  step oslom_real_more "$PY -u real_more_detection_v3.py 10"
  step oslom_real_big  "$PY -u run_oslom_big.py --workers 6 --hours 48 --out results/big_dl/results_big"
  step score_big_oslom "$PY -u score_big.py"
  step score_onmi_oslom "$PY -u score_big.py --onmi-only"
fi

# --- synthetic benchmarks: first round (earlier tables) ----------------------------------------------
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered commonality_family; do
  step run1_$EXP "$PY -u run.py $EXP --workers $W"
done
for EXP in rings_uniform rings_hetero rings_pref; do
  step run1x_$EXP "$PY -u run.py $EXP --workers $W --methods slpa_x,demon_x,angel_x,cpm_x,scan_x"
done
for EXP in lfr lfr_slices lfr_clustered; do
  step run1slpa_$EXP "$PY -u run.py $EXP --workers $W --methods slpa100"
done
for EXP in rings_uniform rings_hetero rings_pref lfr; do
  step selftune_$EXP "$PY -u eval_self_tune.py $EXP 3"
done
echo "$(date '+%F %T') QUEUE END" >> "$RUNS/queue.log"
