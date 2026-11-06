#!/bin/bash
# YouTube: score the 128 GB machine's HyperCommon covers, then the rivals on the same graph.
cd "/d/Python Projects/hypercommon/paper/bench"
PY=../../.venv/Scripts/python.exe
for S in full induced; do
  $PY -u real_detection_v4.py youtube $S > results/real_v4_youtube_$S.log 2>&1
  echo "youtube $S done $(date +%H:%M) exit $?" >> results/queue9_state.txt
done
