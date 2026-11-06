#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
$PY -u real_family_detection.py dblp > results/real_family_dblp2.log 2>&1
echo "family done $(date +%H:%M)" > results/queue5_state.txt
$PY -u runtime.py > results/runtime.log 2>&1
echo "runtime done $(date +%H:%M)" >> results/queue5_state.txt
