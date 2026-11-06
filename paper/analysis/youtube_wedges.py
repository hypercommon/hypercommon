"""
Wedge counts of the SNAP YouTube graph, with and without degree-ratio pruning: since
J(u,v) <= (min(d_u,d_v)+1)/(max(d_u,d_v)+1), an edge whose degree ratio is below the
smallest threshold t_min of interest can never lie in an admitted triple, so the wedges
through it need not be scored. Writes figures/data/youtube_wedges.csv.
"""
import gzip
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
e = np.loadtxt(gzip.open(os.path.join(HERE, "snap", "com-youtube.ungraph.txt.gz"), "rt"), comments="#", dtype=np.int64)
n = e.max() + 1
d = np.bincount(e.ravel(), minlength=n)
rows = [dict(t_min=0.0, edges_kept=1.0, wedges=int((d * (d - 1) // 2).sum()), max_degree=int(d.max()))]
for tmin in (0.02, 0.05, 0.1):
    a, b = d[e[:, 0]], d[e[:, 1]]
    ok = (np.minimum(a, b) + 1) / (np.maximum(a, b) + 1) >= tmin
    dp = np.bincount(e[ok].ravel(), minlength=n)
    rows.append(dict(t_min=tmin, edges_kept=round(float(ok.mean()), 4), wedges=int((dp * (dp - 1) // 2).sum()),
                     max_degree=int(dp.max())))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "..", "figures", "data", "youtube_wedges.csv"), index=False)
print(out.to_string())
