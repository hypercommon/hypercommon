"""
Merge the rows of leiden_cpm_lo (Leiden-CPM at the sixteen resolutions below the old floor
0.002, added after review; queue_review.sh) into leiden_cpm, so that L-CPM's grid is
1e-4 .. 0.5 everywhere the analysis reads it. Renames the method field in place in every
results2/<exp>.csv and <exp>_dense.csv; idempotent (a second run finds nothing to rename).

    python merge_cpm_lo.py [--dry]
"""
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
dry = "--dry" in sys.argv
for path in sorted(glob.glob(os.path.join(HERE, "results2", "*.csv"))):
    name = os.path.basename(path)
    if any(name.endswith(s) for s in ("_summary2.csv", "_paired2.csv", "_selected2.csv",
                                      "_frozen_params2.csv", "_jctrl.csv")):
        continue
    with open(path, newline="") as fh:
        head = fh.readline()
        if "method" not in head.split(","):
            continue
        body = fh.read()
    n = body.count(",leiden_cpm_lo,")
    if not n:
        continue
    col = head.rstrip("\r\n").split(",").index("method")
    rows = list(csv.reader([head] + body.splitlines(keepends=False)))
    k = 0
    for r in rows[1:]:
        if len(r) > col and r[col] == "leiden_cpm_lo":
            r[col] = "leiden_cpm"
            k += 1
    print(f"{name}: {k} rows renamed")
    if dry:
        continue
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        csv.writer(fh, lineterminator="\r\n" if head.endswith("\r\n") else "\n").writerows(rows)
    os.replace(tmp, path)
