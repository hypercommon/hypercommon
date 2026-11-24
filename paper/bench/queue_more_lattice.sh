#!/usr/bin/env bash
# BigClam and ego-splitting on the lattice-like real networks (real_lattice_detection.py; rows
# appended to results/real_lattice_<kind>.csv; resumable).   W=2 bash queue_more_lattice.sh
PY=${PY:-/d/tmp/hcenv314/Scripts/python}
W=${W:-2}
cd "$(dirname "$0")"
for K in protein sketch text tssb spatial transit streets tracks image geo hic; do
  $PY real_lattice_detection.py $K $W bigclam,egosplit >> results/real_lattice_more.log 2>&1
  echo "$K done $(date)" >> results/real_lattice_more_queue.log
done
