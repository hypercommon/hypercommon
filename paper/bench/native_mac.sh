#!/bin/bash
# The steps of queue_mac.sh that its first Mac run missed (native defaults; DEMON on the
# swapped geometry), run beside it with the same step names. Usage, from paper/bench:
#   PY=~/hc/mamba/envs/hc/bin/python WORKERS=4 bash native_mac.sh
cd "$(dirname "$0")"
PY=${PY:-python}
W=${WORKERS:-4}
RUNS=${RUNS:-$HOME/hc/runs}
step() {
  local name=$1; shift
  [ -e "$RUNS/done/$name" ] && { echo "skip $name"; return; }
  echo "$(date '+%F %T') start $name" >> "$RUNS/queue.log"
  if (eval "$@") >> "$RUNS/logs/$name.log" 2>&1; then
    touch "$RUNS/done/$name"; echo "$(date '+%F %T') done  $name" >> "$RUNS/queue.log"
  else
    echo "$(date '+%F %T') FAIL  $name (exit $?)" >> "$RUNS/queue.log"
  fi
}
NATIVE=leiden_g1,leiden_cpm_dens,linkcomm_pd
for EXP in rings_hetero_overlap rings_uniform rings_hetero rings_pref lfr lfr_slices lfr_clustered rings_geometry nulls; do
  step run2native_$EXP "$PY -u run2.py $EXP --workers $W --methods $NATIVE"
done
step run2_rings_geometry_swap_demon "$PY -u run2.py rings_geometry_swap --workers $W --methods demon"
echo "$(date '+%F %T') NATIVE END" >> "$RUNS/queue.log"
