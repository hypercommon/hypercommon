"""
Wedge counts W = sum_v C(d_v, 2) of real networks without ground truth, against m ln m,
with and without the degree-ratio pruning of Sec. "Cost" at t_min = 0.02 (as
youtube_wedges.py): HyperCommon's cost is O(n + m + W). Every graph is read as simple and
undirected (directions, multiple edges and self-loops dropped). SNAP edge lists in
analysis/snap/more/ (https://snap.stanford.edu/data/), plus the three SNAP graphs of the
main comparison in analysis/snap/.

    python wedge_scaling.py [--table-only]
Writes figures/data/wedge_scaling.csv and the Supplemental table
sections/table_wedges.tex (--table-only: the table from the CSV, no SNAP files needed).
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = os.path.join(HERE, "snap")
NETS = [  # (file, kind)
    ("ca-GrQc.txt.gz", "collaboration"), ("ca-HepTh.txt.gz", "collaboration"),
    ("ca-HepPh.txt.gz", "collaboration"), ("ca-AstroPh.txt.gz", "collaboration"),
    ("ca-CondMat.txt.gz", "collaboration"),
    ("email-Enron.txt.gz", "communication"), ("email-EuAll.txt.gz", "communication"),
    ("wiki-Talk.txt.gz", "communication"),
    ("facebook_combined.txt.gz", "social"), ("twitter_combined.txt.gz", "social"),
    ("wiki-Vote.txt.gz", "social"), ("soc-Epinions1.txt.gz", "social"),
    ("soc-Slashdot0902.txt.gz", "social"), ("loc-brightkite_edges.txt.gz", "social"),
    ("loc-gowalla_edges.txt.gz", "social"), ("com-lj.ungraph.txt.gz", "social"),
    ("as-skitter.txt.gz", "internet"), ("as-caida20071105.txt.gz", "internet"),
    ("oregon1_010526.txt.gz", "internet"), ("p2p-Gnutella31.txt.gz", "internet"),
    ("web-Google.txt.gz", "web"), ("web-Stanford.txt.gz", "web"),
    ("web-NotreDame.txt.gz", "web"), ("web-BerkStan.txt.gz", "web"),
    ("roadNet-CA.txt.gz", "road"), ("roadNet-PA.txt.gz", "road"),
    ("cit-HepPh.txt.gz", "citation"), ("cit-Patents.txt.gz", "citation"),
    ("amazon0601.txt.gz", "co-purchase"),
]
MAIN = [("com-amazon.ungraph.txt.gz", "co-purchase"), ("com-dblp.ungraph.txt.gz", "collaboration"),
        ("com-youtube.ungraph.txt.gz", "social")]
T_MIN = 0.02
CSV = os.path.join(HERE, "..", "figures", "data", "wedge_scaling.csv")
TEX = os.path.join(HERE, "..", "sections", "table_wedges.tex")
SHOW = {"facebook_combined": "ego-Facebook", "twitter_combined": "ego-Twitter",
        "loc-brightkite_edges": "loc-Brightkite", "loc-gowalla_edges": "loc-Gowalla",
        "com-lj": "com-LiveJournal", "oregon1_010526": "oregon1-010526",
        "as-caida20071105": "as-CAIDA", "soc-Slashdot0902": "soc-Slashdot"}


def edges(path):
    e = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=[0, 1],
                    dtype=np.int64, engine="c").to_numpy()
    e = e[e[:, 0] != e[:, 1]]
    e.sort(axis=1)
    _, inv = np.unique(e.ravel(), return_inverse=True)       # relabel to 0..n-1
    e = inv.reshape(-1, 2).astype(np.int64)
    n = int(e.max()) + 1
    key = np.unique(e[:, 0] * n + e[:, 1])
    return np.c_[key // n, key % n], n


def row(name, kind, path):
    e, n = edges(path)
    m = len(e)
    d = np.bincount(e.ravel(), minlength=n)
    W = int((d * (d - 1) // 2).sum())
    a, b = d[e[:, 0]], d[e[:, 1]]
    ok = (np.minimum(a, b) + 1) / (np.maximum(a, b) + 1) >= T_MIN
    dp = np.bincount(e[ok].ravel(), minlength=n)
    Wp = int((dp * (dp - 1) // 2).sum())
    mlm = m * np.log(m)
    return dict(network=name, kind=kind, n=n, m=m, max_degree=int(d.max()),
                mean_degree=round(2 * m / n, 2), W=W, W_over_mlnm=round(W / mlm, 2),
                W_pruned=Wp, W_pruned_over_mlnm=round(Wp / mlm, 2))


def table(df):
    """The Supplemental table: the 29 networks, then the three of the main comparison."""
    main = {f.split(".")[0] for f, _ in MAIN}
    lines = [r"\begin{table}[h]",
             r"\caption{Wedges of thirty-two SNAP networks, read as simple undirected graphs: "
             r"nodes $n$, edges $m$, largest degree $d_{\max}$, and the number of wedges "
             r"$W=\sum_v\binom{d_v}{2}$ against $m\ln m$, without and with the degree-ratio "
             r"pruning of Sec.~\ref{sec:cost} at $t_{\min}=0.02$, which is exact for every "
             r"threshold $t\ge0.02$. The last three are the graphs of Table~\ref{tab:real}; "
             r"the others have no ground truth.}",
             r"\label{tab:S-wedges}",
             r"\begin{ruledtabular}",
             r"\begin{tabular}{llrrrrr}",
             r"Network & Kind & $n$ & $m$ & $d_{\max}$ & $W/(m\ln m)$ & Pruned \\",
             r"\hline"]
    order = ([r for r in df.itertuples() if r.network not in main]
             + [r for r in df.itertuples() if r.network in main])
    for i, r in enumerate(order):
        if r.network in main and order[i - 1].network not in main:
            lines.append(r"\hline")
        num = [f"{x:,}".replace(",", r"\,") for x in (r.n, r.m, r.max_degree)]
        lines.append(f"{SHOW.get(r.network, r.network)} & {r.kind} & {' & '.join(num)} & "
                     f"{r.W_over_mlnm:.2f} & {r.W_pruned_over_mlnm:.2f} " + "\\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    with open(TEX, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    if "--table-only" in sys.argv:
        table(pd.read_csv(CSV))
        return
    rows = []
    for f, kind in NETS + MAIN:
        path = os.path.join(SNAP, "more", f) if (f, kind) in NETS else os.path.join(SNAP, f)
        if not os.path.exists(path):
            print("missing", f)
            continue
        rows.append(row(f.split(".")[0], kind, path))
        print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(CSV, index=False)
    print(out.to_string(index=False))
    table(out)


if __name__ == "__main__":
    main()
