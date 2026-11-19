# Sliding-window graphs from the Time Series Segmentation Benchmark (TSSB)

120 graph folders: 24 TSSB series, each as one window-space k-NN graph
(`<series>_knn_k8`), two time-banded graphs (`<series>_band_k8`, `<series>_band_k16`) and
two time-only controls (`<series>_band_k8_null`, `<series>_band_k16_null`). Nodes are
sliding windows; the ground-truth groups are the annotated segments. Built by
`paper/analysis/tssb_build.py`. No community detection was run while building these graphs
and no rule was changed after looking at any result; only the structural numbers below were
computed.

## Sources

| What | URL |
|---|---|
| Repository | https://github.com/ermshaua/time-series-segmentation-benchmark |
| Pinned commit | `e0c7a1f6cee6874e60cbacd3b8601f2de259606e` (2025-03-31, head of `main` on 2026-11-13) |
| Series list, window sizes, change points | https://raw.githubusercontent.com/ermshaua/time-series-segmentation-benchmark/e0c7a1f6cee6874e60cbacd3b8601f2de259606e/tssb/datasets/desc.txt |
| Series values (one file per series) | `.../e0c7a1f6cee6874e60cbacd3b8601f2de259606e/tssb/datasets/<name>.txt` (same prefix) |
| Licence | `.../e0c7a1f6cee6874e60cbacd3b8601f2de259606e/LICENSE` |

- Benchmark: A. Ermshaus, P. Schäfer, U. Leser, "ClaSP: parameter-free time series
  segmentation", *Data Mining and Knowledge Discovery* 37 (2023),
  doi:10.1007/s10618-023-00923-x; and P. Schäfer, A. Ermshaus, U. Leser, "ClaSP - Time
  Series Segmentation", CIKM 2021. 75 univariate series, each built by concatenating
  z-normalised UCR/UEA classification series grouped by class; the change points are the
  concatenation offsets, and each series has a hand-chosen window size.
- `desc.txt` (2,164 bytes, SHA-256 `edcfa27e821c3a83...`): one line per series,
  `name,window_size,cp_1,...,cp_m`. Series files: one float per line.
- Licence: BSD 3-Clause, copyright (c) 2021 Arik Ermshaus (the repository's `LICENSE`).
  The underlying series come from the UCR/UEA time series classification archive
  (http://timeseriesclassification.com/); cite the TSSB papers above and the archive.
- Access: the plain-text files were downloaded from `raw.githubusercontent.com` at the
  pinned commit (no login, no account). The benchmark's own Python package
  (`tssb/utils.py`) was not imported or run; the files are parsed by `tssb_build.py`.

## Selection

1. Of the 75 series, 6 have a single segment (no change point) and so only one group;
   they are dropped. Rule: keep series with 2-9 segments. This leaves 69 series.
2. Six folders per series for all 69 would give 414 folders, far over the ~120 budget.
   So only k = 8 is used for the window-space graph (graph A), and a seeded subset of
   24 series is taken: `random.Random(2026).sample(sorted(names), 24)`, kept in sorted
   order. 24 x 5 = 120 folders. The rule was fixed before any graph was built;
   `selection.json` lists the 69 eligible and the 24 chosen names.
3. All chosen series fit the size limits: n = 1,175-2,996 windows, m <= 27,332 edges.

Chosen (number of segments): Adiac (4), CricketZ (5), DiatomSizeReduction (3), Fish (7),
GunPoint (2), InlineSkate (3), Lightning2 (2), Lightning7 (5), Mallat (4), Meat (3),
MedicalImages (2), MelbournePedestrian (7), OSULeaf (3), Plane (7),
ProximalPhalanxOutlineAgeGroup (3), ProximalPhalanxOutlineCorrect (4), ShapesAll (5),
SonyAIBORobotSurface1 (2), Symbols (4), ToeSegmentation1 (2), UWaveGestureLibraryAll (5),
UWaveGestureLibraryX (5), UWaveGestureLibraryY (5), Yoga (2).

## Construction

Let the series be x[0..len-1], w the benchmark window size (10-50 here) and
c_1 < ... < c_m the change points.

1. **Nodes.** Stride s = max(1, ceil(len / 3000)) (1-6 here). Window i covers samples
   [i*s, i*s + w), for i = 0..n-1 with n = floor((len - w)/s) + 1. Ids are in time order.
2. **Normalisation.** Each window is z-normalised on its own: subtract its mean, divide by
   its population standard deviation (ddof = 0). A window with std <= 1e-8 would become
   the zero vector (none occurred). This is the usual choice for shape-based window
   comparison (matrix profile, ClaSP) and removes level and scale offsets.
3. **Distance.** Squared Euclidean distance between z-normalised windows, rounded to 9
   decimals so that floating-point near-ties are exact ties.
4. **Graph A, `knn`.** Each window i keeps the k = 8 windows j != i with the smallest
   distance (ties: smaller time gap |i - j|, then smaller j). The undirected, unweighted
   edge set is the union over i (so degree >= k). No exclusion zone: overlapping
   neighbours in time are allowed but are rarely chosen (see `frac_windows_overlap`).
   This graph ignores time order.
5. **Graph B, `band`.** The Hi-C recipe: candidates for i are only the windows j with
   0 < |i - j| <= W, with W = max(ceil(5w / s), 2k) windows, that is about five window
   lengths either side in time, but never fewer than 2k candidates per side (the floor
   applies to Fish and Yoga at k = 8, and to the 11 series with 5w/s < 32 at k = 16).
   Each window keeps the k = 8 or 16 candidates with
   the smallest distance (same tie rule); union over i.
6. **Control, `<band graph>_null`.** Same nodes and groups; each window keeps its k
   nearest windows in time (|i - j|, ties to the smaller id), union over i. This is a pure
   band with no signal, built as in `hic_null.py`. One control per series and k; the
   `knn_k8` graph is compared with `band_k8_null`.
7. **Groups.** Segment t covers samples [c_{t-1}, c_t), with c_0 = 0 and c_{m+1} = len.
   Window i is in segment t if it shares at least one sample with it
   (i*s < c_t and i*s + w > c_{t-1}). A window that straddles a change point is therefore
   in both segments: about (w - 1)/s windows per change point (1-65 nodes per series in
   more than one group). Groups with fewer than 3 nodes, and exact duplicate node sets,
   are dropped (none were). Every node is in at least one group.

## Files

Each folder has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v` with u < v, ids 0..n-1, no header,
  no self-loops, no duplicates.
- `groups.txt` - one segment per line (in time order), space-separated sorted node ids.
- `nodes.tsv` - header `id start end segments`; window = samples [start, end) of the
  series; `segments` = comma-separated indices (0-based) of the segments it touches.

`groups.txt` and `nodes.tsv` are identical in all five folders of a series.
`summary.json` has one row per non-control graph (72 rows) with, among others: `n`, `m`,
`mean_degree`, `clustering` (networkx `average_clustering`), `n_components`, `n_groups`,
`group_size_min/median/max`, `frac_nodes_in_group`, `nodes_in_gt1_group`,
`frac_edges_within_group` (edges whose ends share at least one group),
`random_pair_frac_within_group` (the same fraction over all node pairs),
`median_gap_windows` (median |i - j| over edges), `frac_gap_le_k`, `frac_windows_overlap`
(edges joining windows that share samples), and the same graph numbers for the matching
control with prefix `null_`. `selection.json` records the selection.

## Reproduce

Python 3.12.10, numpy 2.5.4, networkx 3.7 (own venv). From `paper/`, with the download
kept in its own folder outside the repository:

    python -I analysis/tssb_build.py fetch --raw <raw_dir>
    python -I analysis/tssb_build.py build --raw <raw_dir> --out analysis/data_tssb
    python -I analysis/tssb_build.py table --out analysis/data_tssb

`fetch` gets 79 small text files (3.4 MB) at the pinned commit; `build` takes about three
minutes and is deterministic (two runs gave byte-identical output).

## Counts

"within" = fraction of edges whose two ends share at least one group; "random pairs" =
the same fraction over all node pairs. Medians over the 24 series:

| graph type | mean degree | clustering | within | within, control | random pairs | frac edges with gap <= k | graphs with > 1 component |
|---|---|---|---|---|---|---|---|
| knn k=8 | 11.2 | 0.418 | 0.538 | 1.000 | 0.308 | 0.008 | 5 of 24 |
| band k=8 | 10.3 | 0.449 | 0.982 | 1.000 | 0.308 | 0.216 | 1 of 24 |
| band k=16 | 19.8 | 0.505 | 0.976 | 1.000 | 0.308 | 0.374 | 0 of 24 |

Controls: mean degree 8.0 (k = 8) and 16.0 (k = 16), clustering 0.64-0.65 and 0.70-0.71,
always connected.

Per graph (w = window, s = stride, W = band half-width in windows, comps = connected
components, sizes = group size min-median-max):

| graph | n | w | s | W | m | mean degree | clustering | comps | groups | sizes | within | within, control | clustering, control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Adiac_knn_k8 | 1399 | 10 | 1 | - | 7220 | 10.32 | 0.665 | 7 | 4 | 176-339-572 | 0.770 | 1.000 | 0.646 |
| Adiac_band_k8 | 1399 | 10 | 1 | 50 | 7611 | 10.88 | 0.478 | 1 | 4 | 176-339-572 | 0.959 | 1.000 | 0.646 |
| Adiac_band_k16 | 1399 | 10 | 1 | 50 | 14047 | 20.08 | 0.525 | 1 | 4 | 176-339-572 | 0.956 | 1.000 | 0.705 |
| CricketZ_knn_k8 | 2708 | 10 | 1 | - | 15143 | 11.18 | 0.347 | 1 | 5 | 477-543-609 | 0.401 | 1.000 | 0.644 |
| CricketZ_band_k8 | 2708 | 10 | 1 | 50 | 14333 | 10.59 | 0.324 | 1 | 5 | 477-543-609 | 0.970 | 1.000 | 0.644 |
| CricketZ_band_k16 | 2708 | 10 | 1 | 50 | 26867 | 19.84 | 0.392 | 1 | 5 | 477-543-609 | 0.973 | 1.000 | 0.702 |
| DiatomSizeReduction_knn_k8 | 1274 | 20 | 1 | - | 6346 | 9.96 | 0.594 | 1 | 3 | 345-450-517 | 0.566 | 1.000 | 0.646 |
| DiatomSizeReduction_band_k8 | 1274 | 20 | 1 | 100 | 6207 | 9.74 | 0.466 | 1 | 3 | 345-450-517 | 0.946 | 1.000 | 0.646 |
| DiatomSizeReduction_band_k16 | 1274 | 20 | 1 | 100 | 12007 | 18.85 | 0.518 | 1 | 3 | 345-450-517 | 0.943 | 1.000 | 0.705 |
| Fish_knn_k8 | 2529 | 10 | 4 | - | 13842 | 10.95 | 0.477 | 1 | 7 | 306-364-407 | 0.234 | 0.998 | 0.645 |
| Fish_band_k8 | 2529 | 10 | 4 | 16 | 13723 | 10.85 | 0.510 | 1 | 7 | 306-364-407 | 0.982 | 0.998 | 0.645 |
| Fish_band_k16 | 2529 | 10 | 4 | 32 | 26663 | 21.09 | 0.560 | 1 | 7 | 306-364-407 | 0.965 | 0.994 | 0.703 |
| GunPoint_knn_k8 | 1866 | 10 | 1 | - | 10535 | 11.29 | 0.479 | 1 | 2 | 900-937.5-975 | 0.638 | 1.000 | 0.645 |
| GunPoint_band_k8 | 1866 | 10 | 1 | 50 | 9830 | 10.54 | 0.469 | 1 | 2 | 900-937.5-975 | 0.990 | 1.000 | 0.645 |
| GunPoint_band_k16 | 1866 | 10 | 1 | 50 | 18705 | 20.05 | 0.542 | 1 | 2 | 900-937.5-975 | 0.991 | 1.000 | 0.703 |
| InlineSkate_knn_k8 | 2406 | 10 | 2 | - | 14658 | 12.18 | 0.394 | 1 | 3 | 529-827-1058 | 0.362 | 1.000 | 0.645 |
| InlineSkate_band_k8 | 2406 | 10 | 2 | 25 | 13556 | 11.27 | 0.549 | 1 | 3 | 529-827-1058 | 0.993 | 1.000 | 0.645 |
| InlineSkate_band_k16 | 2406 | 10 | 2 | 32 | 25490 | 21.19 | 0.609 | 1 | 3 | 529-827-1058 | 0.989 | 0.999 | 0.703 |
| Lightning2_knn_k8 | 1175 | 20 | 1 | - | 6725 | 11.45 | 0.387 | 1 | 2 | 398-597-796 | 0.630 | 1.000 | 0.646 |
| Lightning2_band_k8 | 1175 | 20 | 1 | 100 | 6095 | 10.37 | 0.363 | 1 | 2 | 398-597-796 | 0.968 | 1.000 | 0.646 |
| Lightning2_band_k16 | 1175 | 20 | 1 | 100 | 11389 | 19.39 | 0.422 | 1 | 2 | 398-597-796 | 0.971 | 1.000 | 0.706 |
| Lightning7_knn_k8 | 2164 | 50 | 3 | - | 13468 | 12.45 | 0.360 | 1 | 5 | 266-441-654 | 0.323 | 1.000 | 0.645 |
| Lightning7_band_k8 | 2164 | 50 | 3 | 84 | 12722 | 11.76 | 0.337 | 1 | 5 | 266-441-654 | 0.949 | 1.000 | 0.645 |
| Lightning7_band_k16 | 2164 | 50 | 3 | 84 | 25134 | 23.23 | 0.402 | 1 | 5 | 266-441-654 | 0.949 | 1.000 | 0.703 |
| Mallat_knn_k8 | 1980 | 10 | 2 | - | 11482 | 11.60 | 0.581 | 1 | 4 | 384-452-704 | 0.386 | 1.000 | 0.645 |
| Mallat_band_k8 | 1980 | 10 | 2 | 25 | 10944 | 11.05 | 0.481 | 1 | 4 | 384-452-704 | 0.979 | 1.000 | 0.645 |
| Mallat_band_k16 | 1980 | 10 | 2 | 32 | 20440 | 20.65 | 0.542 | 1 | 4 | 384-452-704 | 0.979 | 0.998 | 0.703 |
| Meat_knn_k8 | 1676 | 10 | 2 | - | 8658 | 10.33 | 0.681 | 24 | 3 | 560-560-564 | 0.762 | 1.000 | 0.645 |
| Meat_band_k8 | 1676 | 10 | 2 | 25 | 8136 | 9.71 | 0.433 | 1 | 3 | 560-560-564 | 0.990 | 1.000 | 0.645 |
| Meat_band_k16 | 1676 | 10 | 2 | 32 | 16617 | 19.83 | 0.525 | 1 | 3 | 560-560-564 | 0.984 | 0.999 | 0.704 |
| MedicalImages_knn_k8 | 2465 | 10 | 1 | - | 14647 | 11.88 | 0.354 | 1 | 2 | 742-1237-1732 | 0.649 | 1.000 | 0.645 |
| MedicalImages_band_k8 | 2465 | 10 | 1 | 50 | 13924 | 11.30 | 0.372 | 1 | 2 | 742-1237-1732 | 0.990 | 1.000 | 0.645 |
| MedicalImages_band_k16 | 2465 | 10 | 1 | 50 | 27332 | 22.18 | 0.433 | 1 | 2 | 742-1237-1732 | 0.990 | 1.000 | 0.703 |
| MelbournePedestrian_knn_k8 | 2444 | 10 | 2 | - | 13559 | 11.10 | 0.526 | 1 | 7 | 288-292-436 | 0.496 | 1.000 | 0.645 |
| MelbournePedestrian_band_k8 | 2444 | 10 | 2 | 25 | 10580 | 8.66 | 0.560 | 1 | 7 | 288-292-436 | 0.969 | 1.000 | 0.645 |
| MelbournePedestrian_band_k16 | 2444 | 10 | 2 | 32 | 22678 | 18.56 | 0.520 | 1 | 7 | 288-292-436 | 0.964 | 0.997 | 0.703 |
| OSULeaf_knn_k8 | 2071 | 10 | 1 | - | 11216 | 10.83 | 0.307 | 1 | 3 | 400-782-907 | 0.640 | 1.000 | 0.645 |
| OSULeaf_band_k8 | 2071 | 10 | 1 | 50 | 10247 | 9.90 | 0.295 | 1 | 3 | 400-782-907 | 0.988 | 1.000 | 0.645 |
| OSULeaf_band_k16 | 2071 | 10 | 1 | 50 | 19073 | 18.42 | 0.356 | 1 | 3 | 400-782-907 | 0.986 | 1.000 | 0.703 |
| Plane_knn_k8 | 1886 | 10 | 2 | - | 9941 | 10.54 | 0.632 | 2 | 7 | 166-270-360 | 0.728 | 1.000 | 0.645 |
| Plane_band_k8 | 1886 | 10 | 2 | 25 | 9027 | 9.57 | 0.421 | 1 | 7 | 166-270-360 | 0.970 | 1.000 | 0.645 |
| Plane_band_k16 | 1886 | 10 | 2 | 32 | 17267 | 18.31 | 0.443 | 1 | 7 | 166-270-360 | 0.961 | 0.996 | 0.703 |
| ProximalPhalanxOutlineAgeGroup_knn_k8 | 1996 | 10 | 2 | - | 11066 | 11.09 | 0.448 | 5 | 3 | 360-699-945 | 0.763 | 1.000 | 0.645 |
| ProximalPhalanxOutlineAgeGroup_band_k8 | 1996 | 10 | 2 | 25 | 9235 | 9.25 | 0.624 | 1 | 3 | 360-699-945 | 0.990 | 1.000 | 0.645 |
| ProximalPhalanxOutlineAgeGroup_band_k16 | 1996 | 10 | 2 | 32 | 18022 | 18.06 | 0.504 | 1 | 3 | 360-699-945 | 0.985 | 0.999 | 0.703 |
| ProximalPhalanxOutlineCorrect_knn_k8 | 2996 | 10 | 2 | - | 16894 | 11.28 | 0.411 | 2 | 4 | 388-699-1222 | 0.631 | 1.000 | 0.644 |
| ProximalPhalanxOutlineCorrect_band_k8 | 2996 | 10 | 2 | 25 | 13867 | 9.26 | 0.627 | 1 | 4 | 388-699-1222 | 0.990 | 1.000 | 0.644 |
| ProximalPhalanxOutlineCorrect_band_k16 | 2996 | 10 | 2 | 32 | 27193 | 18.15 | 0.507 | 1 | 4 | 388-699-1222 | 0.987 | 0.999 | 0.702 |
| ShapesAll_knn_k8 | 2131 | 10 | 3 | - | 12461 | 11.69 | 0.457 | 1 | 5 | 427-429-430 | 0.386 | 1.000 | 0.645 |
| ShapesAll_band_k8 | 2131 | 10 | 3 | 17 | 11682 | 10.96 | 0.583 | 1 | 5 | 427-429-430 | 0.988 | 1.000 | 0.645 |
| ShapesAll_band_k16 | 2131 | 10 | 3 | 32 | 22776 | 21.38 | 0.581 | 1 | 5 | 427-429-430 | 0.974 | 0.996 | 0.703 |
| SonyAIBORobotSurface1_knn_k8 | 1381 | 20 | 1 | - | 8031 | 11.63 | 0.399 | 1 | 2 | 420-700-980 | 0.719 | 1.000 | 0.646 |
| SonyAIBORobotSurface1_band_k8 | 1381 | 20 | 1 | 100 | 6992 | 10.13 | 0.386 | 1 | 2 | 420-700-980 | 0.980 | 1.000 | 0.646 |
| SonyAIBORobotSurface1_band_k16 | 1381 | 20 | 1 | 100 | 13284 | 19.24 | 0.439 | 1 | 2 | 420-700-980 | 0.976 | 1.000 | 0.705 |
| Symbols_knn_k8 | 1880 | 10 | 1 | - | 10491 | 11.16 | 0.476 | 1 | 4 | 298-406.5-796 | 0.511 | 1.000 | 0.645 |
| Symbols_band_k8 | 1880 | 10 | 1 | 50 | 10168 | 10.82 | 0.521 | 4 | 4 | 298-406.5-796 | 0.974 | 1.000 | 0.645 |
| Symbols_band_k16 | 1880 | 10 | 1 | 50 | 19912 | 21.18 | 0.570 | 1 | 4 | 298-406.5-796 | 0.970 | 1.000 | 0.703 |
| ToeSegmentation1_knn_k8 | 2761 | 10 | 1 | - | 15897 | 11.52 | 0.364 | 1 | 2 | 1385-1385-1385 | 0.643 | 1.000 | 0.644 |
| ToeSegmentation1_band_k8 | 2761 | 10 | 1 | 50 | 14444 | 10.46 | 0.416 | 1 | 2 | 1385-1385-1385 | 0.995 | 1.000 | 0.644 |
| ToeSegmentation1_band_k16 | 2761 | 10 | 1 | 50 | 27265 | 19.75 | 0.476 | 1 | 2 | 1385-1385-1385 | 0.994 | 1.000 | 0.702 |
| UWaveGestureLibraryAll_knn_k8 | 2817 | 10 | 3 | - | 15229 | 10.81 | 0.315 | 1 | 5 | 525-544-625 | 0.487 | 1.000 | 0.644 |
| UWaveGestureLibraryAll_band_k8 | 2817 | 10 | 3 | 17 | 12919 | 9.17 | 0.380 | 1 | 5 | 525-544-625 | 0.991 | 1.000 | 0.644 |
| UWaveGestureLibraryAll_band_k16 | 2817 | 10 | 3 | 32 | 25238 | 17.92 | 0.409 | 1 | 5 | 525-544-625 | 0.980 | 0.997 | 0.702 |
| UWaveGestureLibraryX_knn_k8 | 2809 | 10 | 1 | - | 15845 | 11.28 | 0.287 | 1 | 5 | 530-550-625 | 0.451 | 1.000 | 0.644 |
| UWaveGestureLibraryX_band_k8 | 2809 | 10 | 1 | 50 | 14259 | 10.15 | 0.347 | 1 | 5 | 530-550-625 | 0.981 | 1.000 | 0.644 |
| UWaveGestureLibraryX_band_k16 | 2809 | 10 | 1 | 50 | 26077 | 18.57 | 0.412 | 1 | 5 | 530-550-625 | 0.975 | 1.000 | 0.702 |
| UWaveGestureLibraryY_knn_k8 | 2809 | 10 | 1 | - | 16188 | 11.53 | 0.292 | 1 | 5 | 530-550-625 | 0.419 | 1.000 | 0.644 |
| UWaveGestureLibraryY_band_k8 | 2809 | 10 | 1 | 50 | 14291 | 10.18 | 0.385 | 1 | 5 | 530-550-625 | 0.979 | 1.000 | 0.644 |
| UWaveGestureLibraryY_band_k16 | 2809 | 10 | 1 | 50 | 26275 | 18.71 | 0.450 | 1 | 5 | 530-550-625 | 0.974 | 1.000 | 0.702 |
| Yoga_knn_k8 | 2661 | 10 | 6 | - | 15424 | 11.59 | 0.424 | 1 | 2 | 1216-1331-1446 | 0.511 | 0.999 | 0.644 |
| Yoga_band_k8 | 2661 | 10 | 6 | 16 | 13095 | 9.84 | 0.494 | 1 | 2 | 1216-1331-1446 | 0.997 | 0.999 | 0.644 |
| Yoga_band_k16 | 2661 | 10 | 6 | 32 | 26291 | 19.76 | 0.534 | 1 | 2 | 1216-1331-1446 | 0.994 | 0.999 | 0.702 |

## Caveats

- **The groups are few and huge.** 2-7 segments per series, of 166-1,732 windows each,
  covering every node. They are a near-partition into long time intervals, not many small
  communities. A random node pair already shares a group 15-60% of the time (median 31%).
- **"within" is not a signal test here.** Any graph whose edges are short in time puts
  almost every edge inside a segment: the time-only control scores 0.994-1.000. The band
  graphs (0.94-1.00) score slightly *below* their controls in all 48 cases, because they choose some
  longer edges by shape and a few of those cross a change point. So, unlike Hi-C, the
  control cannot show a gain in `frac_edges_within_group`; it only shows how much of a
  detection score could come from time contiguity alone.
- **The knn graphs mix segments.** Only 23-77% of knn edges (median 54%) lie within a
  segment, against 31% for random pairs. The segments in a TSSB series come from classes
  of one UCR dataset and often share local shapes at w = 10-50 samples, so many nearest
  windows sit in another segment. Five knn graphs are disconnected (Meat 24 components,
  Adiac 7, ProximalPhalanxOutlineAgeGroup 5, Plane 2, ProximalPhalanxOutlineCorrect 2);
  of the band graphs only Symbols_band_k8 is (4 components).
- **Little time-local structure inside the band.** Only a median 22% (k = 8) and 37%
  (k = 16) of band edges join windows at most k positions apart, so the band graphs are
  not close to a 1D lattice. Most series are quasi-periodic, and the nearest window in
  shape is often one period away, not the next window. Clustering of the band graphs
  (0.29-0.63) is below that of the time-only lattice (0.64 at k = 8, 0.70 at k = 16).
- Windows that straddle a change point belong to both segments (1-65 nodes per series),
  so the overlaps are small and only at segment borders.
- The window sizes are the benchmark's hand-chosen values (10-50 in the chosen series),
  not tuned; the stride and band rules depend only on the series length and w.
- UWaveGestureLibraryX and UWaveGestureLibraryY have the same length, window size and
  change points (so identical groups and controls); their values differ.
- The 24 series are a seeded random subset of 69; no rule was set using detection
  results.
