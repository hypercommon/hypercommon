#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
while powershell -NoProfile -Command "if (@(Get-CimInstance Win32_Process -Filter \"name='python.exe'\" | Where-Object { \$_.CommandLine -match 'real_more_detection_v2' }).Count -gt 0) { exit 0 } else { exit 1 }"; do sleep 60; done
$PY -u real_detection_v3.py amazon full > results/real3_amazon_full.log 2>&1
$PY -u real_detection_v3.py dblp full > results/real3_dblp_full.log 2>&1
$PY -u real_detection_v3.py amazon induced > results/real3_amazon_induced.log 2>&1
$PY -u real_detection_v3.py dblp induced > results/real3_dblp_induced.log 2>&1
$PY -u runtime.py > results/runtime.log 2>&1
echo "queue6 done $(date +%H:%M)" > results/queue6_state.txt
