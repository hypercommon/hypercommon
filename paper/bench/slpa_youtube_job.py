"""SLPA on the YouTube subgraph induced by the ground-truth nodes, the one big_job.py
task left undone there (its task_simple failed earlier, in clique percolation, with a
MemoryError above 100 GB). Writes results/big_dl/results_big/youtube_induced_slpa.*"""
import os
import sys
import time

import networkx as nx
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import big_job as B                   # noqa: E402
from slpa import slpa_multi           # noqa: E402

B.DATA = os.path.join(HERE, "..", "analysis", "snap")
out = os.path.join(HERE, "results", "big_dl", "results_big")
e, ids = B.load("youtube", "induced")
n = len(ids)
A = B.adjacency(e, n)
G = nx.Graph()
G.add_nodes_from(range(n))
G.add_edges_from(map(tuple, e.tolist()))
st = B.Store(out, "youtube", "induced", "slpa")
np.random.seed(1)
t1 = time.time()
for r, cov in slpa_multi(G, 100, [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45]):
    st.add(r, *B.to_arrays(cov), A, ids, time.time() - t1)
    print(r, len(cov), round(time.time() - t1, 1), flush=True)
    t1 = time.time()
st.save()
print("done")
