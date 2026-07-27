"""Smoke test for pcrit_tbest_omega — n=400 config."""
from experiments.ring_lattice import pcrit_tbest_omega as exp

exp.N         = 400
exp.Z         = 16
exp.RING_SIZE = 100
exp.P_STEP    = 0.01  # full 101 p values
exp.T_STEP    = 0.01  # full 101 t values
exp.N_RUNS    = 3

if __name__ == "__main__":
    exp.run_experiment(out_dir="results/_smoke_pcrit_tbest_omega", seed=7, n_runs=exp.N_RUNS)
