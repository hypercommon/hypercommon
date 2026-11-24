#!/usr/bin/env bash
# After review: (1) Leiden-CPM at sixteen resolutions below its grid (leiden_cpm_lo, 1e-4 .. 0.0017)
# on every benchmark experiment where L-CPM runs; (2) the same in the robustness sweep (step
# 0.01) for the gradient families, the only ones with rings large enough (600 nodes, k = 16) to be
# split by CPM at gamma = 0.002 (the old floor); merge_cpm_lo.py then renames the rows to
# leiden_cpm; (3) the control leiden_j (Leiden-RB on Jaccard edge weights) on the ring and LFR
# benchmarks, rows to results2/<exp>_jctrl.csv. Resumable.
#   bash queue_review.sh
PY=${PY:-/d/tmp/hcenv314/Scripts/python}
W=${W:-20}
cd "$(dirname "$0")"
LOG=results2/review_queue.log
for EXP in rings_hetero rings_hetero_overlap rings_uniform rings_pref lfr lfr_slices lfr_clustered \
           rings_geometry nulls rings_isolated rings_hubs rings_split rings_geometry_swap rgg; do
  $PY run2.py $EXP --methods leiden_cpm_lo --workers $W >> results2/${EXP}_cpmlo.log 2>&1
  echo "$EXP cpm_lo done $(date)" >> $LOG
done
DENSE_P=1 $PY run2.py rings_hetero --methods leiden_cpm_lo --families gradient10 --workers $W \
  >> results2/rings_hetero_dense_cpmlo.log 2>&1
DENSE_P=1 $PY run2.py rings_hetero_overlap --methods leiden_cpm_lo --families gradient10_o00,gradient10_o10 \
  --workers $W >> results2/rings_hetero_overlap_dense_cpmlo.log 2>&1
echo "dense gradient cpm_lo done $(date)" >> $LOG
# ego-splitting at step 0.01 (not 0.05) where it is the runner-up in p_0.6 (unequal size, o_n = 0.10)
DENSE_P=1 DENSE_STEP=0.01 $PY run2.py rings_hetero_overlap --methods egosplit --families sizes_only_k8_o10 \
  --workers $W >> results2/rings_hetero_overlap_dense_ego01.log 2>&1
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered; do
  RUN2_SUFFIX=_jctrl $PY run2.py $EXP --methods leiden_j --workers $W >> results2/${EXP}_jctrl.log 2>&1
  echo "$EXP leiden_j done $(date)" >> $LOG
done
echo "ALL DONE $(date)" >> $LOG
