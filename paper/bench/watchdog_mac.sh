#!/bin/bash
# Memory watchdog for the runs on the Mac (macOS). Every 5 s it reads each of our Python
# processes' memory as top reports it (resident plus compressed), and force-kills
#   - any of them above PER_GB, and
#   - the largest of them once the machine's swap passes SWAP_GB.
# A killed multiprocessing worker would leave its pool waiting forever, so the pool's
# parent and its other workers are killed with it: that step fails and the queue goes on.
# Only processes of the env in ~/hc/mamba are touched, and the OSLOM2 binary in
# ~/tools/OSLOM2, which is killed alone: its Python caller records the failed run and goes
# on. Log: ~/hc/runs/logs/watchdog.log
PER_GB=${PER_GB:-30}
SWAP_GB=${SWAP_GB:-12}
LOG=$HOME/hc/runs/logs/watchdog.log
ENV=/hc/mamba/envs/

ours() { ps -o command= -p "$1" 2>/dev/null | awk '{print $1}' | grep -q -e "${ENV}.*/bin/python" -e "/tools/OSLOM2/oslom_undir$"; }
oslom() { ps -o command= -p "$1" 2>/dev/null | awk '{print $1}' | grep -q "/tools/OSLOM2/oslom_undir$"; }

kill_tree() {   # kill_tree <pid> <reason>
  local pid=$1 parent
  parent=$(ps -o ppid= -p "$pid" | tr -d ' ')
  if oslom "$pid"; then
    echo "$(date '+%F %T') KILL oslom $pid ($2)" >> "$LOG"
    kill -9 "$pid"
  elif [ -n "$parent" ] && ours "$parent"; then
    echo "$(date '+%F %T') KILL pool $parent and workers (worker $pid: $2)" >> "$LOG"
    pkill -9 -P "$parent"; kill -9 "$parent"
  else
    echo "$(date '+%F %T') KILL $pid ($2): $(ps -o command= -p "$pid" | cut -c1-120)" >> "$LOG"
    pkill -9 -P "$pid"; kill -9 "$pid"
  fi
}

to_gb() {       # top's MEM column (e.g. 100G, 512M, 900K) in GB
  awk -v m="$1" 'BEGIN{u=substr(m,length(m)); v=m+0; if(u=="G")print v; else if(u=="M")print v/1024; else if(u=="K")print v/1048576; else print v/1073741824}'
}

echo "$(date '+%F %T') watchdog start: per process ${PER_GB} GB, swap ${SWAP_GB} GB" >> "$LOG"
while true; do
  big_pid=; big_gb=0
  while read -r pid mem cmd; do
    case "$cmd" in python*|oslom_undir*) ;; *) continue ;; esac
    ours "$pid" || continue
    gb=$(to_gb "${mem%[+-]}")
    if awk -v a="$gb" -v b="$PER_GB" 'BEGIN{exit !(a>b)}'; then
      kill_tree "$pid" "${gb} GB > ${PER_GB} GB"; continue
    fi
    if awk -v a="$gb" -v b="$big_gb" 'BEGIN{exit !(a>b)}'; then big_pid=$pid; big_gb=$gb; fi
  done < <(top -l 1 -o mem -n 40 -stats pid,mem,command 2>/dev/null | awk 'NR>1 && $1 ~ /^[0-9]+$/')
  swap=$(sysctl -n vm.swapusage | awk '{for(i=1;i<=NF;i++) if($i=="used") {v=$(i+2); sub("M","",v); print v/1024}}')
  if [ -n "$big_pid" ] && awk -v a="$swap" -v b="$SWAP_GB" 'BEGIN{exit !(a>b)}'; then
    kill_tree "$big_pid" "swap ${swap} GB > ${SWAP_GB} GB, largest at ${big_gb} GB"
  fi
  sleep 5
done
