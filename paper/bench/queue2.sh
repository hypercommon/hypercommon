#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
running() { powershell -NoProfile -Command "@(Get-CimInstance Win32_Process -Filter \"name='python.exe'\" | Where-Object { \$_.CommandLine -match 'run.py commonality_family' }).Count"; }
while [ "$(running | tr -d '\r')" != "0" ]; do sleep 60; done
$PY -u run.py rings_hetero --workers 20 > results_rings_hetero2.log 2>&1
$PY -u run.py lfr_clustered --workers 20 > results_lfr_clustered.log 2>&1
$PY -u run.py rings_pref --workers 20 > results_rings_pref.log 2>&1
$PY -u run.py commonality_family --workers 20 > results_family2.log 2>&1
$PY -u runtime.py > results/runtime.log 2>&1
$PY -u run.py rings_uniform --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_uniform.log 2>&1
$PY -u run.py rings_hetero --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_hetero.log 2>&1
$PY -u run.py rings_pref --workers 20 --methods slpa_x,demon_x,angel_x,cpm_x,scan_x > results_ext_pref.log 2>&1
