#!/bin/bash
# Main runs, third pass: queue10 and queue11 ran at the same time and ran out of
# memory, so the jobs that failed are rerun here one at a time (run2.py resumes).
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
ALL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,demon,angel,scan,scan_j,leiden_g1,leiden_cpm_dens,linkcomm_pd
CTRL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,leiden_g1,leiden_cpm_dens,linkcomm_pd
NAT=leiden_g1,leiden_cpm_dens,linkcomm_pd
W=16
$PY -u run2.py rings_hetero_overlap --workers $W --methods $ALL > results2/rings_hetero_overlap_v3c.log 2>&1
echo "q12 hetero_overlap done $(date +%F_%H:%M)" >> results2/queue10_state.txt
for EXP in rings_pref lfr lfr_clustered; do
  $PY -u run2.py $EXP --workers $W --methods $NAT > results2/${EXP}_native_c.log 2>&1
  echo "q12 $EXP native done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
for EXP in rings_hubs rings_split rings_geometry_swap; do
  $PY -u run2.py $EXP --workers $W --methods $CTRL > results2/${EXP}_c.log 2>&1
  echo "q12 $EXP done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
echo "queue12 done $(date +%F_%H:%M)" >> results2/queue10_state.txt
