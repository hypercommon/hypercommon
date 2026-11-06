#!/bin/bash
# Mop-up after queue12: rerun the tasks that ran out of memory when two queues ran at
# once (run2.py now counts a MemoryError row as not done).
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
ALL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,demon,angel,scan,scan_j,leiden_g1,leiden_cpm_dens,linkcomm_pd
CTRL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
until grep -q "queue12 done" results2/queue10_state.txt; do sleep 60; done
for EXP in rings_hetero rings_uniform rings_pref lfr lfr_slices; do
  $PY -u run2.py $EXP --workers 16 --methods $ALL > results2/${EXP}_mop.log 2>&1
done
$PY -u run2.py rgg --workers 16 --methods $CTRL > results2/rgg_mop.log 2>&1
echo "queue13 done $(date +%F_%H:%M)" >> results2/queue10_state.txt
