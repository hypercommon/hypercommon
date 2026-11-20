#!/bin/sh
cd "$(dirname "$0")"
for k in more hic sketch image tracks geo protein transit streets tssb spatial text; do
  ${PYTHON:-python} real_tree_detection.py $k 18 || echo FAILED $k
done
echo ALL DONE
