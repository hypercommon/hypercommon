"""BigClam port (bench/bigclam.py) against SNAP's own binary (examples/bigclam, one thread),
on the smaller networks at several K: best-match F1 and Omega against the recorded groups.
SNAP built under WSL in ~/tools/snap (examples/bigclam, make with a conda g++); its random
seed is fixed (TAGMFast(G, 10, 10)), the port runs seeds 0-2.

    python check_bigclam_snap.py football,polbooks 5,12,30 > ../figures/data/bigclam_vs_snap.csv"""
import os, subprocess, sys
import numpy as np

B = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, B); sys.path.insert(1, os.path.join(B, "..", "analysis"))
os.chdir(B)
import cover_metrics as M
from real_more import load
from bigclam import bigclam_covers

W = r"D:\tmp\bigclam_check"            # seen from WSL as /mnt/d/tmp/bigclam_check (below)
os.makedirs(W, exist_ok=True)
NETS = sys.argv[1].split(",")
KS = [int(k) for k in sys.argv[2].split(",")]
SEEDS = range(3)


def score(cov, truth, universe, nodes_u):
    r = [set(c) & universe for c in cov]
    r = [c for c in r if c]
    return M.f1_best_match(r, truth, universe), M.omega(r, truth, nodes_u)


print("network,K,impl,seed,f1,omega,n_comms", flush=True)
for net in NETS:
    G, truth = load(net)
    universe = set().union(*truth); nodes_u = list(universe)
    lab = {u: i for i, u in enumerate(G.nodes())}
    inv = {i: u for u, i in lab.items()}
    ef = os.path.join(W, f"{net}.txt")
    with open(ef, "w", newline="\n") as fh:
        for a, b in G.edges():
            if a != b:
                fh.write(f"{lab[a]}\t{lab[b]}\n")
    for K in KS:
        if K > G.number_of_nodes():
            continue
        pre = f"/mnt/d/tmp/bigclam_check/{net}_K{K}_"
        cmd = (f"cd ~/tools/snap/examples/bigclam && ./bigclam -i:/mnt/d/tmp/bigclam_check/{net}.txt "
               f"-c:{K} -nt:1 -o:{pre} > /dev/null 2>&1")
        subprocess.run(["wsl", "-e", "bash", "-lc", cmd], check=True)
        cov = []
        with open(os.path.join(W, f"{net}_K{K}_cmtyvv.txt")) as fh:
            for line in fh:
                ids = line.split()
                if ids:
                    cov.append({inv[int(x)] for x in ids})
        f1, om = score(cov, truth, universe, nodes_u)
        print(f"{net},{K},snap,10,{f1:.4f},{om:.4f},{len(cov)}", flush=True)
        for s in SEEDS:
            (_, c), = list(bigclam_covers(G, [K], s))
            f1, om = score(c, truth, universe, nodes_u)
            print(f"{net},{K},port,{s},{f1:.4f},{om:.4f},{len(c)}", flush=True)
