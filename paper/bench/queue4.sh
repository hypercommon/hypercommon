#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
$PY -u run.py lfr --workers 12 --methods slpa100 > results_slpa100_lfr.log 2>&1
$PY -u run.py lfr_slices --workers 12 --methods slpa100 > results_slpa100_slices.log 2>&1
$PY -u run.py lfr_clustered --workers 12 --methods slpa100 > results_slpa100_clustered.log 2>&1
$PY -u run.py rings_hetero_overlap --workers 16 > results_hetero_overlap.log 2>&1
$PY -u real_family_detection.py dblp > results/real_family_dblp2.log 2>&1
$PY -u runtime.py > results/runtime.log 2>&1
echo queue4 finished
