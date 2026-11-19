#!/usr/bin/env bash
# Robustness sweep, OSLOM and the stochastic block model (WSL): p = 0, 0.01, ..., 1 minus the
# main grid at step DENSE_STEP (0.05: adds p = 0.85 ... 1), test graphs only. Rows to results2/<exp>_dense_{oslom,sbm}.csv;
# resumable (a rerun skips finished tasks). W = workers per method (two methods at once).
#   bash queue_dense_wsl.sh [--dry]
PY=${PY:-~/tools/mamba/envs/gt/bin/python}
export OSLOM_BIN=${OSLOM_BIN:-$HOME/tools/OSLOM2/oslom_undir}
export DENSE_P=1
export DENSE_STEP=${DENSE_STEP:-0.05}   # OSLOM and SBM: about 100 s per graph, so the 0.05 grid only
W=${W:-3}
cd "$(dirname "$0")"
for EXP in rings_uniform rings_hetero rings_hetero_overlap rings_pref; do
  if [ "$1" = "--dry" ]; then
    $PY run_oslom.py $EXP --dry
    $PY run_sbm.py $EXP --dry
    continue
  fi
  $PY run_oslom.py $EXP --workers $W >> results2/${EXP}_dense_oslom.log 2>&1 &
  $PY run_sbm.py $EXP --workers $W >> results2/${EXP}_dense_sbm.log 2>&1 &
  wait
  echo "$EXP done $(date)" >> results2/dense_wsl_queue.log
done
