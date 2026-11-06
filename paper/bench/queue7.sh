#!/bin/bash
# Benchmark main runs: every experiment with methods2.py, sequentially.
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
mkdir -p results2
for EXP in rings_isolated rings_hetero rings_hetero_overlap rings_pref rings_uniform lfr lfr_slices lfr_clustered rings_geometry; do
  $PY -u run2.py $EXP --workers 16 > results2/$EXP.log 2>&1
  echo "$EXP done $(date +%H:%M)" >> results2/queue7_state.txt
done
echo "queue7 done $(date +%H:%M)" >> results2/queue7_state.txt
