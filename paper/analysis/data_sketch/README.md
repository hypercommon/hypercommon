# Point graphs from hand-drawn sketches, with pen strokes as ground truth

40 graphs: 20 drawings from the Google Quick, Draw! dataset (4 each from the
categories bicycle, cat, face, flower, house), each with k = 8 and k = 16,
plus an (x, y)-only control `<name>_null` for every graph (80 folders in all).
Nodes are points resampled at a fixed arc-length step along the pen strokes;
edges join a point to its nearest points in (x, y, scaled time); the groups
are the strokes, so every group is a long, thin chain. Built by
`paper/analysis/sketch_build.py` (run with `python -I`). No community
detection was run, and no rule was changed after looking at any number below;
only the structural statistics below were computed.

## Sources

| What | URL |
|---|---|
| Raw drawings, one file per category | `https://storage.googleapis.com/quickdraw_dataset/full/raw/<category>.ndjson` |
| Dataset page and format | https://github.com/googlecreativelab/quickdraw-dataset (README read at commit `5fe6c0a910b3732bc3db7639d3c9e7c287617f2f`, head of `master` on 2026-11-13) |
| Licence | https://creativecommons.org/licenses/by/4.0/ |

- Data: "The Quick, Draw! Dataset", Google Creative Lab (2017). Drawings made
  by players of the game https://quickdraw.withgoogle.com, stored as strokes
  of `[x, y, t]` lists, with x, y in pixels and t the time in ms since the
  drawing's first point.
- Licence: the dataset README states that the data is made available by
  Google, Inc. under the Creative Commons Attribution 4.0 International
  licence (CC BY 4.0). Credit "Quick, Draw! Dataset, Google" when using these
  graphs. Changes made here: coordinates rescaled, strokes resampled, graphs
  built (see below).
- Download: public Google Cloud Storage bucket, no login, curl, 2026-11-13.
  The full raw files are 438-593 MB, so only the first 20,000,000 bytes of
  each were fetched (`curl -r 0-19999999`); the last, cut-off line of each
  prefix is never read. The files were only parsed with `json.loads`; nothing
  from them was executed.

| category | full file (bytes) | SHA-256 of the 20,000,000-byte prefix |
|---|---|---|
| bicycle | 520,247,189 | `60d92e512a457b6e9de6a4995d08b2c9370c4bc3137a7f8f6cdeca7737d32631` |
| cat | 521,081,210 | `04a97d78682399f3e6552799b06b6f2fc894d6ef232fe69ab40a2405dc023c58` |
| face | 593,069,574 | `b631e0bce85325829e6c6ee66cb1abae0292e9712b75fe6e5966303840bd719f` |
| flower | 553,084,310 | `817d829470daeabe6c8e58464782a67398cbb522af817e9abbbaaf8599e84ce2` |
| house | 437,886,394 | `b447f7c3173c17e3c23d092f6093a1591eb3be27d2a0a15b31d1e9ccf97e3a6f` |

## Selection (fixed before building)

1. Categories: bicycle, cat, face, flower, house (chosen by hand, before
   looking at the data, as everyday objects that are usually drawn with
   several separate strokes).
2. From each prefix file, the first 4000 complete lines (drawings).
3. Eligible: `recognized == true`; at least 6 strokes; `t` non-decreasing
   inside every stroke and from each stroke's end to the next stroke's start;
   median pen speed > 0; 300 <= n <= 1500 nodes after resampling.
4. Eligible drawings sorted by `key_id` (as integers). One
   `rng = random.Random(2026)`; for the categories in the order above,
   `rng.sample(eligible, 4)`, then sorted by `key_id` again. No drawing was
   swapped or re-drawn.

Counts among the 4000 scanned drawings (also in `selection.json`; a rule is
checked only if the earlier ones passed, in the order of step 3):

| category | eligible | unrecognized | < 6 strokes | bad timing | n out of range |
|---|---|---|---|---|---|
| bicycle | 2826 | 120 | 1052 | 0 | 2 |
| cat | 2991 | 665 | 340 | 3 | 1 |
| face | 1933 | 324 | 1722 | 2 | 19 |
| flower | 529 | 149 | 3319 | 0 | 3 |
| house | 574 | 96 | 3327 | 2 | 1 |

## Construction

1. **Normalise.** Per drawing, shift x and y so their minimum is 0, and scale
   both by the same factor so that max(x range, y range) = 255 units (as in
   the official "simplified" files). `t` stays in ms.
2. **Resample.** Per stroke, consecutive points with the same (x, y) are
   merged (the first is kept). The polyline is resampled at arc length
   s = 0, 2, 4, ..., 2*floor(L/2) (fixed step 2.0 units, L = stroke length);
   x, y and t are interpolated linearly in s between raw points. A stroke of
   length 0 gives one node. Node ids run in drawing order (stroke by stroke,
   along each stroke), so ids also run in time order. x, y are rounded to 4
   decimals and t to 3, and all edges are computed from these rounded values.
3. **Time scaling (exact).** v = median, over all consecutive raw point pairs
   inside a stroke with dt > 0 (coordinates after step 1, before the merging
   of step 2), of |dxy| / dt, in units per ms. Then tau = v * t. So 1 ms
   becomes the distance the pen covers in 1 ms at its median speed: along a
   stroke drawn at median speed, tau grows exactly as fast as arc length
   (one time unit = one spatial unit). v is per drawing, stored as
   `median_speed_units_per_ms` in `summary.json` (median 0.28, range
   0.17-0.83).
4. **Real graph `<name>_k8`, `<name>_k16` (exact rule).** Each node i keeps
   the k nodes j != i with the smallest Euclidean distance in (x, y, tau);
   squared distances rounded to 9 decimals, ties to smaller j. The undirected,
   unweighted edge set is the union over i, so every degree is >= k. No
   self-loops; each edge written once as `u<TAB>v` with u < v.
5. **Control `<name>_k<k>_null`.** Same nodes, same `groups.txt` and
   `nodes.tsv`; same rule with the distance in (x, y) only (no time).
6. **Groups.** The nodes of each stroke with >= 3 nodes, one line per stroke
   in drawing order. Groups are disjoint. Nodes of shorter strokes (stroke
   length < 4 units) are in no group: 8 nodes in all, in 5 drawings.

`frac_within` = fraction of edges whose two ends are in the same group.

## Reproduce

Needs numpy and networkx (we reused the image-graph venv read-only: Python
3.12.10, numpy 2.5.4, networkx 3.7). Fetch the five prefixes into their own
empty folder:

    curl -r 0-19999999 -o <cat>.part.ndjson https://storage.googleapis.com/quickdraw_dataset/full/raw/<cat>.ndjson

then from `paper/`:

    python -I analysis/sketch_build.py build --raw <folder> --out analysis/data_sketch
    python -I analysis/sketch_build.py table --out analysis/data_sketch   # table below

`build` takes about 2 minutes on one core; two runs gave byte-identical files.

## Files

Each folder `<category><key_id>_k<K>/` and `<category><key_id>_k<K>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v.
- `groups.txt` - one stroke per line, space-separated sorted node ids, size >= 3.
- `nodes.tsv` - header `id x y t`; x, y in normalised units (0-255), t in ms
  since the drawing's first point (the real graph uses tau = v * t, v from
  `summary.json`). No stroke id.

`summary.json` has one row per real graph, with the control's numbers under
`null_*`: n, m, mean/min/max degree, clustering (networkx
`average_clustering`), n_components, `frac_within`, `frac_edges_same_stroke`
(also counts strokes with < 3 nodes), `frac_edges_consecutive` (edges between
neighbouring samples of one stroke), n_strokes, n_strokes_small,
n_nodes_in_small_strokes, n_groups, group size min/median/mean/max,
`frac_nodes_in_group`, stroke length min/median/max (units, groups only),
total ink length, duration, median stroke duration, median and min pen-up gap
(ms), v, and the drawing's `key_id`, `countrycode` and `timestamp`.
`selection.json` has the scan counts and the chosen key_ids.

## Counts

20 drawings, n = 428-847 (median 632.5; 12,506 nodes in all), 6-17 strokes
(median 8.5), 4-17 groups (median 8). 176 groups in all: sizes 4-426 nodes
(median 49.5, mean 71.0); mean group size per drawing 40-141. Stroke length
(groups only): per-drawing median 23.6-238.0 units (median of medians 96.3),
longest stroke 145-851 units; total ink 849-1689 units per drawing. Drawing
time 4.0-16.8 s; median pen-up gap 138-878 ms, i.e. 54-254 tau units.

Median over the 20 drawings (min-max in brackets):

| quantity | k = 8 | k = 8 null | k = 16 | k = 16 null |
|---|---|---|---|---|
| m | 2766 (1932-3486) | 2684 (1849-3626) | 5620 (3915-7167) | 5531 (3842-7382) |
| mean degree | 8.66 (8.23-9.03) | 8.55 (8.16-8.97) | 17.72 (16.67-18.29) | 17.51 (16.46-18.51) |
| clustering | 0.712 (0.672-0.739) | 0.655 (0.627-0.704) | 0.787 (0.727-0.828) | 0.698 (0.667-0.735) |
| components | 7.5 (3-14) | 2 (1-7) | 6 (3-12) | 1 (1-3) |
| frac_within | 0.997 (0.978-1.000) | 0.907 (0.776-0.997) | 0.987 (0.943-1.000) | 0.833 (0.688-0.993) |
| consecutive edges | 0.228 (0.217-0.241) | 0.231 (0.218-0.242) | 0.112 (0.107-0.119) | 0.113 (0.106-0.120) |

Per graph ("within" = frac_within; "consec." = fraction of edges between
neighbouring samples of one stroke; stroke length in units):

| graph | n | strokes (groups) | group size min-median-max | stroke length median | m | m null | mean deg | mean deg null | within | within null | consec. | consec. null |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bicycle4890906597523456_k8 | 461 | 9 (8) | 8-33-148 | 65.2 | 2035 | 1922 | 8.83 | 8.34 | 0.992 | 0.931 | 0.222 | 0.235 |
| bicycle4890906597523456_k16 | 461 | 9 (8) | 8-33-148 | 65.2 | 4120 | 3972 | 17.87 | 17.23 | 0.946 | 0.858 | 0.110 | 0.114 |
| bicycle5527324114026496_k8 | 719 | 11 (11) | 7-44-186 | 87.3 | 3067 | 3075 | 8.53 | 8.55 | 0.995 | 0.866 | 0.231 | 0.230 |
| bicycle5527324114026496_k16 | 719 | 11 (11) | 7-44-186 | 87.3 | 6265 | 6283 | 17.43 | 17.48 | 0.979 | 0.834 | 0.113 | 0.113 |
| bicycle6150134419161088_k8 | 788 | 13 (13) | 15-41-151 | 81.4 | 3460 | 3534 | 8.78 | 8.97 | 1.000 | 0.776 | 0.224 | 0.219 |
| bicycle6150134419161088_k16 | 788 | 13 (13) | 15-41-151 | 81.4 | 7073 | 7293 | 17.95 | 18.51 | 0.996 | 0.688 | 0.110 | 0.106 |
| bicycle6178682009812992_k8 | 481 | 8 (8) | 12-43-155 | 85.0 | 2090 | 2056 | 8.69 | 8.55 | 1.000 | 0.909 | 0.226 | 0.230 |
| bicycle6178682009812992_k16 | 481 | 8 (8) | 12-43-155 | 85.0 | 4224 | 4177 | 17.56 | 17.37 | 0.979 | 0.846 | 0.112 | 0.113 |
| cat4802956002066432_k8 | 630 | 9 (9) | 19-37-314 | 72.1 | 2813 | 2757 | 8.93 | 8.75 | 1.000 | 0.887 | 0.221 | 0.225 |
| cat4802956002066432_k16 | 630 | 9 (9) | 19-37-314 | 72.1 | 5695 | 5766 | 18.08 | 18.30 | 1.000 | 0.804 | 0.109 | 0.108 |
| cat5201136883597312_k8 | 611 | 9 (9) | 9-40-335 | 78.5 | 2720 | 2678 | 8.90 | 8.77 | 1.000 | 0.907 | 0.221 | 0.225 |
| cat5201136883597312_k16 | 611 | 9 (9) | 9-40-335 | 78.5 | 5498 | 5613 | 18.00 | 18.37 | 0.987 | 0.830 | 0.109 | 0.107 |
| cat6339492900241408_k8 | 680 | 17 (17) | 5-23-163 | 44.2 | 2912 | 3043 | 8.56 | 8.95 | 0.979 | 0.804 | 0.228 | 0.218 |
| cat6339492900241408_k16 | 680 | 17 (17) | 5-23-163 | 44.2 | 5900 | 6233 | 17.35 | 18.33 | 0.957 | 0.699 | 0.112 | 0.106 |
| cat6489727114936320_k8 | 573 | 13 (13) | 13-47-73 | 93.1 | 2531 | 2462 | 8.83 | 8.59 | 1.000 | 0.977 | 0.221 | 0.227 |
| cat6489727114936320_k16 | 573 | 13 (13) | 13-47-73 | 93.1 | 5083 | 5134 | 17.74 | 17.92 | 0.990 | 0.869 | 0.110 | 0.109 |
| face5488597970976768_k8 | 657 | 6 (5) | 42-60-426 | 118.6 | 2897 | 2690 | 8.82 | 8.19 | 0.997 | 0.997 | 0.225 | 0.242 |
| face5488597970976768_k16 | 657 | 6 (5) | 42-60-426 | 118.6 | 5847 | 5449 | 17.80 | 16.59 | 0.993 | 0.993 | 0.111 | 0.119 |
| face5791035743535104_k8 | 847 | 6 (6) | 9-33.5-402 | 66.4 | 3486 | 3626 | 8.23 | 8.56 | 1.000 | 0.812 | 0.241 | 0.232 |
| face5791035743535104_k16 | 847 | 6 (6) | 9-33.5-402 | 66.4 | 7061 | 7382 | 16.67 | 17.43 | 0.987 | 0.744 | 0.119 | 0.114 |
| face5948591019393024_k8 | 498 | 7 (4) | 30-56.5-351 | 112.5 | 2078 | 2045 | 8.35 | 8.21 | 0.986 | 0.984 | 0.236 | 0.240 |
| face5948591019393024_k16 | 498 | 7 (4) | 30-56.5-351 | 112.5 | 4195 | 4168 | 16.85 | 16.74 | 0.985 | 0.978 | 0.117 | 0.118 |
| face6349802449141760_k8 | 635 | 10 (10) | 6-50.5-262 | 99.4 | 2708 | 2625 | 8.53 | 8.27 | 0.993 | 0.987 | 0.231 | 0.238 |
| face6349802449141760_k16 | 635 | 10 (10) | 6-50.5-262 | 99.4 | 5545 | 5313 | 17.46 | 16.73 | 0.985 | 0.958 | 0.113 | 0.118 |
| flower4517225065086976_k8 | 428 | 8 (8) | 4-57.5-87 | 113.8 | 1932 | 1849 | 9.03 | 8.64 | 0.990 | 0.837 | 0.217 | 0.227 |
| flower4517225065086976_k16 | 428 | 8 (8) | 4-57.5-87 | 113.8 | 3915 | 3842 | 18.29 | 17.95 | 0.986 | 0.754 | 0.107 | 0.109 |
| flower5498212439818240_k8 | 473 | 8 (8) | 6-61.5-93 | 121.8 | 2033 | 2022 | 8.60 | 8.55 | 0.991 | 0.890 | 0.229 | 0.230 |
| flower5498212439818240_k16 | 473 | 8 (8) | 6-61.5-93 | 121.8 | 4158 | 4154 | 17.58 | 17.56 | 0.984 | 0.813 | 0.112 | 0.112 |
| flower6494518385836032_k8 | 735 | 6 (6) | 112-119.5-144 | 238.0 | 3135 | 3130 | 8.53 | 8.52 | 1.000 | 0.906 | 0.233 | 0.233 |
| flower6494518385836032_k16 | 735 | 6 (6) | 112-119.5-144 | 238.0 | 6476 | 6487 | 17.62 | 17.65 | 1.000 | 0.846 | 0.113 | 0.112 |
| flower6625890043691008_k8 | 750 | 14 (14) | 4-12.5-248 | 23.6 | 3306 | 3211 | 8.82 | 8.56 | 0.978 | 0.912 | 0.223 | 0.229 |
| flower6625890043691008_k16 | 750 | 14 (14) | 4-12.5-248 | 23.6 | 6668 | 6644 | 17.78 | 17.72 | 0.943 | 0.832 | 0.110 | 0.111 |
| house4832121245925376_k8 | 675 | 6 (5) | 45-100-221 | 198.3 | 2876 | 2866 | 8.52 | 8.49 | 0.997 | 0.883 | 0.233 | 0.233 |
| house4832121245925376_k16 | 675 | 6 (5) | 45-100-221 | 198.3 | 5923 | 5806 | 17.55 | 17.20 | 0.997 | 0.824 | 0.113 | 0.115 |
| house4933089958559744_k8 | 807 | 11 (11) | 9-60-170 | 119.5 | 3470 | 3436 | 8.60 | 8.52 | 1.000 | 0.899 | 0.229 | 0.232 |
| house4933089958559744_k16 | 807 | 11 (11) | 9-60-170 | 119.5 | 7167 | 7006 | 17.76 | 17.36 | 0.990 | 0.845 | 0.111 | 0.114 |
| house4963337613869056_k8 | 548 | 6 (6) | 36-59-262 | 117.0 | 2366 | 2326 | 8.64 | 8.49 | 1.000 | 0.918 | 0.229 | 0.233 |
| house4963337613869056_k16 | 548 | 6 (6) | 36-59-262 | 117.0 | 4853 | 4805 | 17.71 | 17.54 | 1.000 | 0.802 | 0.112 | 0.113 |
| house5675020489588736_k8 | 510 | 6 (5) | 22-74-237 | 146.2 | 2216 | 2081 | 8.69 | 8.16 | 0.996 | 0.985 | 0.227 | 0.242 |
| house5675020489588736_k16 | 510 | 6 (5) | 22-74-237 | 146.2 | 4520 | 4197 | 17.73 | 16.46 | 0.996 | 0.959 | 0.112 | 0.120 |

## Caveats

- **Many real graphs fall apart into one piece per stroke.** Pen-up gaps
  (median 54-254 tau units per drawing) are far longer than the 2-unit
  sample step, so in (x, y, tau) a point's 8 or 16 nearest points are nearly
  always on its own stroke. The real graphs have 3-14 components at k = 8
  (median 7.5) and 3-12 at k = 16 (median 6). In 9 of the 20 drawings at
  k = 8 (2 of 20 at k = 16: cat4802956002066432, house4963337613869056) the
  components are exactly the strokes (frac_within = 1 and components =
  strokes), so the ground truth is just the connected components there. The
  controls are mostly connected (median 2 components at k = 8, 1 at k = 16)
  and are the harder case.
- Every stroke is a 1D chain with the same local structure because of the
  fixed arc-length step: about 23% of edges at k = 8 (11% at k = 16) join
  neighbouring samples, in real and control graphs alike.
- In some drawings one stroke holds most of the ink (face outline, house
  body): the largest group is up to 70% of n (face5948591019393024; 65% in
  face5488597970976768). There the control's frac_within is about as high as
  the real one's (face5488597970976768, face5948591019393024,
  face6349802449141760, house5675020489588736: null 0.96-0.997 at k = 8 and
  k = 16).
- Strokes shorter than 4 units give < 3 nodes and are not groups, so a few
  drawings with >= 6 strokes have only 4-5 groups (face5948591019393024: 7
  strokes, 4 groups).
- Node ids run in time order, so each group is a contiguous id range. A
  method that uses node order could exploit this; shuffle ids if that
  matters.
- Only the first 4000 drawings of each file were scanned; they are not
  documented to be a uniform sample of the category. The chosen drawings
  date from January-March 2017.
- 80 folders (20 drawings x 2 k x real/control) is more than the 60 of the
  image family, but the graphs are small (n <= 847, m <= 7,382).
