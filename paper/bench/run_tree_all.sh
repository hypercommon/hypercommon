#!/bin/sh
# all ring and LFR settings, one after another
cd "$(dirname "$0")"
for e in rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr lfr_slices lfr_clustered; do
  ${PYTHON:-python} run_tree.py $e --workers 18 || exit 1
done
echo ALL DONE
