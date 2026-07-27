"""Smoke test for overlap_algos_experiment — tiny config, validates whole pipeline.

Uses n=200, z=8, ring_size=50 (4 rings), 2 overlaps, p_step=0.05 (20 steps),
AVG_SEARCH=2, AVG_FINAL=2.

Run:
  ./.venv/Scripts/python.exe -m experiments.ring_lattice._smoke_overlap_algos
"""

from experiments.ring_lattice import overlap_algos_experiment as exp

# Override module constants
exp.N            = 200
exp.Z            = 8
exp.RING_SIZE    = 50
exp.P_STEP       = 0.05
exp.AVG_FINAL    = 2
exp.AVG_SEARCH   = 2
exp.COARSE_T_STEP = 0.1
exp.FINE_T_STEP   = 0.05
exp.FINE_T_WINDOW = 0.05

if __name__ == "__main__":
    exp.run_experiment(
        out_dir="results/_smoke_overlap_algos",
        seed=7,
        overlaps=[0.0, 0.04],
    )
