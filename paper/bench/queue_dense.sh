#!/usr/bin/env bash
# Robustness sweep, all main methods (Windows, Git Bash): p = 0, 0.01, ..., 1 minus the main
# grid, test graphs only (DENSE_P=1, run2.py). Rows to results2/<exp>_dense.csv; resumable
# (a rerun skips finished tasks). Then the planted-structure measures (robust_structure.py).
# Run together with queue_dense_wsl.sh; W + 2 x (its W) should stay at or below 22.
#   bash queue_dense.sh [--dry]
PY=${PY:-/d/tmp/hcenv314/Scripts/python}
export DENSE_P=1
W=${W:-16}
# the main methods at step 0.01; BigClam and ego-splitting run at 0.05 (queue_more.sh), OSLOM and the SBM
# at 0.05 in WSL (queue_dense_wsl.sh)
M=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,demon,angel,scan,scan_j
cd "$(dirname "$0")"
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref; do
  if [ "$1" = "--dry" ]; then
    $PY run2.py $EXP --methods $M --dry
    continue
  fi
  $PY run2.py $EXP --methods $M --workers $W >> results2/${EXP}_dense.log 2>&1
  echo "$EXP done $(date)" >> results2/dense_queue.log
done
[ "$1" = "--dry" ] || $PY robust_structure.py --workers $W >> results2/robust_structure.log 2>&1
