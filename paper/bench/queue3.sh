#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
$PY -u run.py rings_uniform --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_uniform.log 2>&1
$PY -u run.py rings_hetero --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_hetero.log 2>&1
$PY -u run.py rings_pref --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_pref.log 2>&1
$PY -u runtime.py > results/runtime.log 2>&1
echo queue3 finished
$PY -u run.py lfr --workers 20 --methods slpa100 > results_slpa100_lfr.log 2>&1
$PY -u run.py lfr_slices --workers 20 --methods slpa100 > results_slpa100_slices.log 2>&1
$PY -u run.py lfr_clustered --workers 20 --methods slpa100 > results_slpa100_clustered.log 2>&1
echo slpa100 finished
$PY -u run.py rings_hetero_overlap --workers 20 > results_hetero_overlap.log 2>&1
$PY -u real_family_detection.py dblp > results/real_family_dblp2.log 2>&1
