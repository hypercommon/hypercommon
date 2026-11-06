#!/bin/bash
# Amazon / DBLP main runs; waits for the 14-network main runs to finish.
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
while powershell -NoProfile -Command "if (@(Get-CimInstance Win32_Process -Filter \"name='python.exe'\" | Where-Object { \$_.CommandLine -match 'real_more_detection_v3' }).Count -gt 0) { exit 0 } else { exit 1 }"; do sleep 120; done
for NET in amazon dblp; do
  for SETTING in full induced; do
    $PY -u real_detection_v4.py $NET $SETTING > results/real4_${NET}_${SETTING}.log 2>&1
    echo "$NET $SETTING done $(date +%H:%M)" >> results/queue8_state.txt
  done
done
echo "queue8 done $(date +%H:%M)" >> results/queue8_state.txt
