# OSLOM on the benchmark graphs (WSL; OSLOM2 built by build_oslom.sh; PY = python with the bench deps)
PY=${PY:-~/tools/mamba/envs/gt/bin/python}
cd "$(dirname "$0")"
for EXP in lfr rings_uniform rings_hetero rings_hetero_overlap rings_pref lfr_slices lfr_clustered nulls rgg; do
  $PY run_oslom.py $EXP --workers ${W:-20} >> results2/${EXP}_oslom.log 2>&1
  echo "$EXP done $(date)" >> results2/oslom_queue.log
done
