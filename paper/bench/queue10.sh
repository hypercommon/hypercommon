#!/bin/bash
# Main runs: the two omitted unequal families at o_n = 0, 0.10 (with every method),
# the native-default methods on every benchmark, and the new control families.
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
ALL=hc,hc_2test,hc_tri,hc_econn,linkcomm,cpm,leiden,leiden_cpm,infomap,lpa,walktrap,slpa,demon,angel,scan,scan_j,leiden_g1,leiden_cpm_dens,linkcomm_pd
NAT=leiden_g1,leiden_cpm_dens,linkcomm_pd
mkdir -p results2
$PY -u run2.py rings_hetero_overlap --workers 20 --methods $ALL > results2/rings_hetero_overlap_v3.log 2>&1
echo "hetero_overlap done $(date +%F_%H:%M)" >> results2/queue10_state.txt
for EXP in rings_hetero rings_uniform rings_pref lfr lfr_slices lfr_clustered rings_geometry; do
  $PY -u run2.py $EXP --workers 20 --methods $NAT > results2/${EXP}_native.log 2>&1
  echo "$EXP native done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
for EXP in rings_hubs rings_split rings_geometry_swap rgg; do
  $PY -u run2.py $EXP --workers 20 --methods $ALL > results2/$EXP.log 2>&1
  echo "$EXP done $(date +%F_%H:%M)" >> results2/queue10_state.txt
done
echo "queue10 done $(date +%F_%H:%M)" >> results2/queue10_state.txt
