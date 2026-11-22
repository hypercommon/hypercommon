#!/usr/bin/env bash
# BigClam and ego-splitting (bigclam.py, egosplit.py), added in a later round: every benchmark
# experiment (rows appended to results2/<exp>.csv), then the robustness sweep at step 0.05 as
# for OSLOM and the SBM (results2/<exp>_dense.csv). Resumable.
#   bash queue_more.sh
PY=${PY:-/d/tmp/hcenv/Scripts/python}
W=${W:-20}
M=bigclam,egosplit
cd "$(dirname "$0")"
for EXP in rings_isolated lfr lfr_slices lfr_clustered nulls rings_uniform rings_hetero rings_hetero_overlap \
           rings_pref rings_geometry rings_hubs rings_split rings_geometry_swap rgg; do
  $PY run2.py $EXP --methods $M --workers $W >> results2/${EXP}_more.log 2>&1
  echo "$EXP done $(date)" >> results2/more_queue.log
done
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref; do
  DENSE_P=1 DENSE_STEP=0.05 $PY run2.py $EXP --methods $M --workers $W >> results2/${EXP}_dense_more.log 2>&1
  echo "$EXP dense done $(date)" >> results2/more_queue.log
done
