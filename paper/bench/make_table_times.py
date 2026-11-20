"""
Table of running times (Supplemental Material, Sec. Cost): results/time_methods.csv
(time_methods.py, Windows and WSL) and results/retime_big_laptop.csv (retime_big.py with
RETIME_OUT=retime_big_laptop.csv), both on the laptop of paper/README.md, Machine.

    python make_table_times.py
Writes ../sections/table_times.tex and ../figures/data/time_methods.csv.
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
LABEL = {"hc": "HC", "linkcomm": "Link", "cpm": "CPM", "leiden": "Leiden", "leiden_cpm": "L-CPM",
         "infomap": "Infomap", "slpa": "SLPA", "oslom": "OSLOM", "sbm": "SBM", "walktrap": "Walktrap",
         "demon": "DEMON", "angel": "ANGEL", "lpa": "LPA", "scan": "SCAN", "scan_j": "SCAN-J"}
ORDER = list(LABEL)
COLS = [("k16_s100_o05", "uniform"), ("two_worlds", "two worlds"), ("on10", "LFR")]


def th(n):
    return f"{int(n):,}".replace(",", "\\,")


def fmt(x):
    if x < 1:
        return f"{x:.2f}"
    if x < 10:
        return f"{x:.1f}"
    return f"{x:.0f}"


def main():
    d = pd.read_csv(os.path.join(HERE, "results", "time_methods.csv"))
    d = d[d.error.isna()]
    t = d.groupby(["method", "family"]).wall_s.median().unstack()
    k = d.groupby("method").n_params.first()
    size = d.groupby("family")[["n_nodes", "n_edges"]].first()
    b = pd.read_csv(os.path.join(HERE, "results", "retime_big_laptop.csv")).set_index(["network", "method"])
    out = t.copy()
    out["n_params"] = k
    for net in ("amazon", "dblp"):
        out[net] = [b.wall_s.get((net, m), float("nan")) for m in out.index]
    out.loc[ORDER].to_csv(os.path.join(HERE, "..", "figures", "data", "time_methods.csv"))

    hc_eq = {net: b.eq_s.get((net, "hc")) for net in ("amazon", "dblp")}
    cap = ("Running time of each method's whole parameter sweep, in seconds: wall time from the "
           "graph to the covers at every value of its grid (HyperCommon: the scoring of all wedges "
           "and the covers at its 30 thresholds), one thread, one method at a time, median of three "
           "runs. Scoring the covers is not timed. Benchmark graphs (seed 0, $p=0.2$ for the rings, "
           "$\\mu=0.3$ for LFR): uniform rings ($s=100$, $k=16$, $o_n=0.05$; "
           f"$n={th(size.n_nodes['k16_s100_o05'])}$, $m={th(size.n_edges['k16_s100_o05'])}$), two worlds "
           f"($o_n=0.05$; $n={th(size.n_nodes['two_worlds'])}$, $m={th(size.n_edges['two_worlds'])}$) and "
           f"overlapping LFR ($o_n=0.1$; $n={th(size.n_nodes['on10'])}$, $m={th(size.n_edges['on10'])}$). "
           "Amazon and DBLP: the full SNAP graphs of Table~\\ref{tab:real}, one run per grid value; "
           "methods not run on them are marked --; there Leiden sweeps 19 resolutions "
           "(\\texttt{bench/big\\_job.py}). The SBM is one fit, its hierarchy levels read off "
           "as its parameter values. Selecting HyperCommon's threshold by overlapping "
           f"modularity adds {hc_eq['amazon']:.0f}\\,s on Amazon and {hc_eq['dblp']:.0f}\\,s on DBLP. "
           "Intel Core Ultra 9 185H laptop; OSLOM and the SBM under WSL on the same machine.")
    lines = ["\\begin{table}[h]", "\\caption{" + cap + "}", "\\label{tab:S-times}",
             "\\begin{ruledtabular}", "\\begin{tabular}{lrrrrrr}",
             "Method & Grid & " + " & ".join(c for _, c in COLS) + " & Amazon & DBLP \\\\", "\\hline"]
    for m in ORDER:
        if m not in out.index:
            continue
        r = out.loc[m]
        big = [fmt(r[n]) if pd.notna(r[n]) else "--" for n in ("amazon", "dblp")]
        grid = 1 if m == "sbm" else int(r.n_params)
        lines.append(f"{LABEL[m]} & {grid} & " + " & ".join(fmt(r[f]) for f, _ in COLS)
                     + " & " + " & ".join(big) + " \\\\")
    lines += ["\\end{tabular}", "\\end{ruledtabular}", "\\end{table}"]
    tex = "\n".join(lines)
    with open(os.path.join(HERE, "..", "sections", "table_times.tex"), "w") as fh:
        fh.write(tex + "\n")
    print(tex)


if __name__ == "__main__":
    main()
