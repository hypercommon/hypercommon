"""Run overlap_algos_experiment with z=8, ring_size=50 (hypothesised to widen
the hypercommon-vs-Leiden gap due to modularity resolution limit + sparser graph).
"""
from experiments.ring_lattice import overlap_algos_experiment as exp

exp.N         = 2000
exp.Z         = 8
exp.RING_SIZE = 50

if __name__ == "__main__":
    # write to a separate sub-dir so we don't confuse with the z=16 run
    exp.run_experiment(out_dir="results/overlap_algos_z8_r50", seed=42)
