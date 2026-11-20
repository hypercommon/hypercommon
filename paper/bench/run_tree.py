"""
The per-branch threshold (hc_tree.py) on the benchmark test graphs, exploratory.
run2.py's task runner with rows to results_tree/<experiment>.csv, test seeds only, so
the paper's result files are untouched.

    python run_tree.py <experiment> [--workers N]
"""
import os
import sys

import run2

if __name__ == "__main__":
    run2.RESULTS = os.path.join(run2.HERE, "results_tree")
    run2.CAL_SEEDS = []
    sys.argv += ["--methods", "hc_tree"]
    run2.main()
