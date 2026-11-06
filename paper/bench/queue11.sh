#!/bin/bash
# Main runs, continued: wait for the running rings_hetero_overlap job (PID 6592), then
# the native-default methods on every benchmark and the control families. The control
# families use the methods that lead on the main benchmarks (DEMON, ANGEL, SCAN, SCAN-J,
# never the most accurate there and the slowest, are left out).
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
ALL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,demon,angel,scan,scan_j,leiden_g1,leiden_cpm_dens,linkcomm_pd
CTRL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
NAT=leiden_g1,leiden_cpm_dens,linkcomm_pd
while tasklist //FI "PID eq 6592" 2>/dev/null | grep -q 6592; do sleep 60; done
echo "hetero_overlap (first pass) done $(date +%F_%H:%M)" >> results2/queue10_state.txt
$PY -u run2.py rings_hetero_overlap --workers 20 --methods $ALL > results2/rings_hetero_overlap_v3b.log 2>&1
echo "hetero_overlap done $(date +%F_%H:%M)" >> results2/queue10_state.txt
for EXP in rings_hetero rings_uniform rings_pref lfr lfr_slices lfr_clustered rings_geometry; do
  $PY -u run2.py $EXP --workers 20 --methods $NAT > results2/${EXP}_native.log 2>&1
  echo "$EXP native done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
for EXP in rings_hubs rings_split rings_geometry_swap rgg; do
  $PY -u run2.py $EXP --workers 20 --methods $CTRL > results2/$EXP.log 2>&1
  echo "$EXP done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
echo "queue11 done $(date +%F_%H:%M)" >> results2/queue10_state.txt
