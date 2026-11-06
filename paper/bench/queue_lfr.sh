#!/bin/bash
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
$PY - <<'P' > results/lfr_pregen.log 2>&1
import experiments as E
n=0
for exp in ("lfr","lfr_slices"):
    X=E.EXPERIMENTS[exp]
    for f,fam in X["families"].items():
        for mu in X["p_grid"]:
            for s in list(range(10))+[1000,1001,1002]:
                E._lfr(fam,mu,s); n+=1
print("generated",n)
P
while tasklist | grep -q python && pgrep -f "run.py rings_hetero" >/dev/null 2>&1 || [ -n "$(powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"name='python.exe'\" | Where-Object { \$_.CommandLine -match 'run.py rings_hetero' } | Select-Object -First 1 -ExpandProperty ProcessId")" ]; do sleep 30; done
$PY -u run.py lfr --workers 20 > results_lfr.log 2>&1
$PY -u run.py lfr_slices --workers 20 > results_lfr_slices.log 2>&1
$PY -u run.py commonality_family --workers 20 > results_family.log 2>&1
$PY -u run.py rings_pref --workers 20 > results_rings_pref.log 2>&1
$PY -u run.py rings_hetero --workers 20 > results_rings_hetero2.log 2>&1
$PY -u run.py lfr_clustered --workers 20 > results_lfr_clustered.log 2>&1
