# Pixel graphs from natural images with human segmentations as ground truth

30 graphs: 15 images from the BSDS500 test split, each with k = 8 and k = 16,
plus a pure-grid control `<name>_null` for every graph (60 folders in all).
Nodes are 6x6-pixel blocks, edges join a block to its most similar-coloured
neighbours, and the groups are the segments drawn by the human annotators.
Built by `paper/analysis/image_build.py` (run with `python -I`).

## Sources

| What | URL |
|---|---|
| BSDS500 archive (images, ground truth, benchmark code) | https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/BSR/BSR_bsds500.tgz |
| Dataset page | https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html |
| Terms (BSDS300 page, which the BSDS500 page extends) | https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/ |

- Dataset: P. Arbelaez, M. Maire, C. Fowlkes, J. Malik, "Contour Detection
  and Hierarchical Image Segmentation", *IEEE TPAMI* 33(5):898-916 (2011),
  doi:10.1109/TPAMI.2010.161. The images and the train/test segmentations
  come from the BSDS300: D. Martin, C. Fowlkes, D. Tal, J. Malik, "A Database
  of Human Segmented Natural Images and its Application to Evaluating
  Segmentation Algorithms and Measuring Ecological Statistics", *ICCV* 2001.
- The old address `http://www.eecs.berkeley.edu/.../BSR/BSR_bsds500.tgz`
  now redirects to a 404; the same path on `www2.eecs.berkeley.edu` (the
  official Berkeley host) works, so no mirror was needed. Downloaded on
  2026-11-13 with curl, no login.
- Archive: 70,763,455 bytes,
  sha256 `97e49d31764f3912f0c4122707d53062ac9e783ba0f095e447a4d53c1a41af8e`,
  md5 `0b17383ff0b21fe43d2f28e624203169`. Berkeley publishes no checksum;
  as a check, the archive holds the expected 200 / 100 / 200 train / val /
  test images, each with a ground-truth `.mat` file of the same id. Only
  `BSR/BSDS500/data/images/test/*.jpg` and
  `BSR/BSDS500/data/groundTruth/test/*.mat` were extracted and read
  (Pillow for the JPEGs, `scipy.io.loadmat` for the `.mat` files). No file
  from the archive was executed.
- Terms: the Berkeley page says the dataset may be downloaded "for
  non-commercial research and educational purposes", asks that work based on
  it cite the ICCV 2001 paper, and asks that segmentation results on the test
  set be shared with the authors. There is no open licence; the photographs
  are from the Corel collection. Cite Martin et al. 2001 and Arbelaez et al.
  2011 when using these graphs.

## Selection (fixed before building)

The 200 test-split ids, sorted as integers, then
`random.Random(2026).sample(ids, 15)`, then sorted again:

```
43051 49024 101027 120003 163096 201080 206097 223060 226022 228076
246009 253016 268048 346016 388067
```

No image was swapped or re-drawn. 15 images rather than ~25 because two k
values plus a control per graph would otherwise exceed the ~60-folder budget
of the downstream benchmark.

## Construction

All rules were fixed before any graph was built. No community detection was
run, and nothing was tuned; only the structural statistics below were computed.

1. **Grid.** Each image (481x321, or 321x481 portrait) is centre-cropped to a
   multiple of 6 in each direction (width x height 480x318, or 318x480).
   Landscape: the last pixel column is cut, and 1 row above and 2 below.
   Portrait: the last pixel row is cut, and 1 column left and 2 right. The
   crop is cut into 6x6 blocks, one node per block: 4240 nodes (53 rows x 80
   columns, or 80 rows x 53 columns for the three portrait images 163096,
   246009 and 388067). Node id = `row * W + col` (row-major).
2. **Colour.** sRGB to CIELAB per pixel (`skimage.color.rgb2lab`, D65,
   2 degree observer), then the mean L, a, b over the block's 36 pixels.
3. **Edges (exact rule).** For node i, the candidates are all nodes j != i
   with Chebyshev grid distance <= 3 (48 candidates inside the image, fewer
   at the border; 15 at a corner). Keep the k candidates with the smallest
   CIELAB distance `|lab_i - lab_j|` (Delta E 1976); ties broken by smaller
   Euclidean grid distance, then by smaller j. The undirected, unweighted
   edge set is the union of these choices over all i, so every degree is
   >= k (>= 15 at the corners when k = 16). No self-loops, each edge written
   once as `u<TAB>v` with u < v.
4. **Groups.** Each annotator's segmentation is mapped to the grid by
   majority: a block gets the label held by most of its 36 pixels (ties to
   the smaller label). Each label of each annotator gives one node set. The
   sets of all annotators (5 or 6 per image) are pooled into one cover:
   exact duplicate node sets are kept once (first occurrence, in the
   annotator order of the `.mat` file), sets with fewer than 3 nodes are
   dropped, and a set holding all n nodes would be dropped (none did). The
   result is an overlapping, nested cover: every node is in one group per
   annotator (minus duplicates), so in 4.9-6.0 groups on average, and coarse
   annotations nest the finer ones.
5. **Control `<name>_null`.** Same nodes, same `groups.txt` and `nodes.tsv`.
   For node i, keep the k nodes j != i with the smallest Euclidean grid
   distance (no radius limit), ties broken by smaller j; union over i. This
   is a pure grid with no colour signal. At k = 8 it is exactly the 8-
   neighbour (king's move) grid (m = 17,099). At k = 16 a node takes its
   12 nodes at distance <= 2 and 4 of the 8 at distance sqrt 5 (the 4 with
   smallest ids), so the union has mean degree 19.8.

Why this rule: it is the Hi-C top-k recipe on a 2D order. The spatial radius
keeps every edge local, so linked nodes share many grid neighbours (the
lattice-like regime), and choosing by colour inside the radius makes edges
avoid object borders, where colour usually changes. The block size was set
to reach about 4000 nodes, and the radius 3 so that k = 16 still leaves a
choice (16 of 48).

## Reproduce

Needs numpy, scipy, Pillow, scikit-image and networkx (we used a separate
venv: Python 3.12.10, numpy 2.5.4, scipy 1.18.1, Pillow 13.0.0,
scikit-image 0.26.0, networkx 3.7). Put the archive in its own empty folder,
extract it there (`tar -xzf BSR_bsds500.tgz`), then from `paper/`:

    python -I analysis/image_build.py build --bsr <folder holding BSR/> --out analysis/data_image
    python -I analysis/image_build.py null  --out analysis/data_image   # optional: rewrite only *_null/

`build` writes all 60 folders and `summary.json` in about a minute. `null`
rebuilds the control edges from `nodes.tsv` and gives byte-identical files.

## Files

Each folder `bsds<ID>_k<K>/` and `bsds<ID>_k<K>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v.
- `groups.txt` - one ground-truth group (one annotator's segment) per line,
  space-separated sorted node ids, size >= 3.
- `nodes.tsv` - header `id row col L a b segment_per_annotator`; `row`,
  `col` are the block's grid position (block covers pixel rows
  `r0 + 6*row .. r0 + 6*row + 5` of the original image, with
  `r0, c0 = crop_offset_rc` in `summary.json`); `L a b` are the mean CIELAB
  values; `segment_per_annotator` is the comma-separated majority label of
  the block under each annotator, in `.mat` order.

`summary.json` has one row per graph, with the control's numbers under
`null_*`. Fields include n, m, mean degree, min/max degree, clustering
(networkx `average_clustering`), n_components, n_groups and the drop counts,
group size min/median/max, `frac_nodes_in_group` (1.0 everywhere),
`groups_per_node_mean`, `group_pairs_jaccard_ge_0_9` (near-duplicate pairs
of groups, mostly the same object drawn by two annotators),
`frac_edges_cheb1` (edges between touching blocks), and two "within" measures:

- `frac_edges_within_group`: fraction of edges whose ends share at least one
  group of the pooled cover (the measure used for the Hi-C graphs);
- `frac_edges_same_segment_annot_mean`: for each annotator, the fraction of
  edges whose ends have the same label, averaged over annotators. With a
  pooled cover the first measure is near 1 for any local graph (one coarse
  annotator is enough), so the second one shows the colour signal better.

## Counts

All graphs have n = 4240 and are connected, and so are all controls. Median
over the 15 images (min-max in brackets):

| quantity | k = 8 | k = 8 null | k = 16 | k = 16 null |
|---|---|---|---|---|
| m | 21,413 (21,065-22,215) | 17,099 | 41,939 (41,371-44,616) | 41,949-41,976 |
| mean degree | 10.10 | 8.07 | 19.78 | 19.79 |
| clustering | 0.496 (0.467-0.532) | 0.442 | 0.544 (0.510-0.580) | 0.551 |
| within (pooled cover) | 0.991 (0.980-0.997) | 0.977 (0.922-0.992) | 0.988 (0.967-0.994) | 0.958 (0.870-0.985) |
| same segment, annotator mean | 0.948 (0.897-0.983) | 0.932 (0.845-0.971) | 0.934 (0.870-0.978) | 0.903 (0.783-0.959) |

Groups per image: 12-203 (median 66), sizes 3 to 3752 nodes.

Per graph ("within" = pooled cover; "same seg." = annotator mean):

| graph | grid | m | mean degree | clustering | groups | size min-median-max | within | within, null | same seg. | same seg., null |
|---|---|---|---|---|---|---|---|---|---|---|
| bsds43051_k8 | 53x80 | 21820 | 10.29 | 0.475 | 12 | 14-666.5-3609 | 0.997 | 0.988 | 0.983 | 0.971 |
| bsds43051_k16 | 53x80 | 42993 | 20.28 | 0.515 | 12 | 14-666.5-3609 | 0.993 | 0.978 | 0.978 | 0.959 |
| bsds49024_k8 | 53x80 | 21350 | 10.07 | 0.502 | 29 | 34-917-1676 | 0.997 | 0.987 | 0.963 | 0.944 |
| bsds49024_k16 | 53x80 | 42108 | 19.86 | 0.544 | 29 | 34-917-1676 | 0.994 | 0.973 | 0.953 | 0.920 |
| bsds101027_k8 | 53x80 | 21368 | 10.08 | 0.475 | 47 | 4-163-2089 | 0.994 | 0.978 | 0.959 | 0.939 |
| bsds101027_k16 | 53x80 | 41674 | 19.66 | 0.523 | 47 | 4-163-2089 | 0.988 | 0.958 | 0.950 | 0.913 |
| bsds120003_k8 | 53x80 | 21494 | 10.14 | 0.532 | 126 | 4-86.5-1048 | 0.991 | 0.963 | 0.936 | 0.874 |
| bsds120003_k16 | 53x80 | 42855 | 20.21 | 0.572 | 126 | 4-86.5-1048 | 0.983 | 0.936 | 0.903 | 0.823 |
| bsds163096_k8 | 80x53 | 21285 | 10.04 | 0.488 | 55 | 3-141-2118 | 0.987 | 0.978 | 0.948 | 0.932 |
| bsds163096_k16 | 80x53 | 41599 | 19.62 | 0.530 | 55 | 3-141-2118 | 0.977 | 0.958 | 0.933 | 0.903 |
| bsds201080_k8 | 53x80 | 21237 | 10.02 | 0.499 | 80 | 3-78.5-2275 | 0.996 | 0.970 | 0.964 | 0.928 |
| bsds201080_k16 | 53x80 | 41564 | 19.61 | 0.549 | 80 | 3-78.5-2275 | 0.993 | 0.951 | 0.953 | 0.899 |
| bsds206097_k8 | 53x80 | 21413 | 10.10 | 0.496 | 87 | 3-96-2440 | 0.995 | 0.992 | 0.935 | 0.924 |
| bsds206097_k16 | 53x80 | 42151 | 19.88 | 0.545 | 87 | 3-96-2440 | 0.992 | 0.985 | 0.918 | 0.893 |
| bsds223060_k8 | 53x80 | 21237 | 10.02 | 0.480 | 72 | 3-77-1481 | 0.980 | 0.975 | 0.954 | 0.941 |
| bsds223060_k16 | 53x80 | 41810 | 19.72 | 0.525 | 72 | 3-77-1481 | 0.975 | 0.956 | 0.942 | 0.913 |
| bsds226022_k8 | 53x80 | 22215 | 10.48 | 0.473 | 81 | 3-36-3209 | 0.989 | 0.989 | 0.906 | 0.920 |
| bsds226022_k16 | 53x80 | 44616 | 21.05 | 0.510 | 81 | 3-36-3209 | 0.982 | 0.978 | 0.890 | 0.891 |
| bsds228076_k8 | 53x80 | 21065 | 9.94 | 0.500 | 66 | 3-39.5-2012 | 0.996 | 0.975 | 0.966 | 0.935 |
| bsds228076_k16 | 53x80 | 41371 | 19.51 | 0.549 | 66 | 3-39.5-2012 | 0.989 | 0.961 | 0.950 | 0.907 |
| bsds246009_k8 | 80x53 | 21229 | 10.01 | 0.512 | 54 | 3-408-1781 | 0.995 | 0.973 | 0.973 | 0.942 |
| bsds246009_k16 | 80x53 | 41939 | 19.78 | 0.550 | 54 | 3-408-1781 | 0.990 | 0.954 | 0.962 | 0.916 |
| bsds253016_k8 | 53x80 | 21542 | 10.16 | 0.509 | 158 | 3-52.5-1786 | 0.991 | 0.963 | 0.897 | 0.857 |
| bsds253016_k16 | 53x80 | 42063 | 19.84 | 0.548 | 158 | 3-52.5-1786 | 0.981 | 0.936 | 0.870 | 0.803 |
| bsds268048_k8 | 53x80 | 21459 | 10.12 | 0.489 | 58 | 3-97.5-3752 | 0.991 | 0.988 | 0.945 | 0.934 |
| bsds268048_k16 | 53x80 | 41919 | 19.77 | 0.528 | 58 | 3-97.5-3752 | 0.989 | 0.979 | 0.934 | 0.907 |
| bsds346016_k8 | 53x80 | 21563 | 10.17 | 0.467 | 60 | 3-115-3067 | 0.986 | 0.977 | 0.937 | 0.930 |
| bsds346016_k16 | 53x80 | 41930 | 19.78 | 0.514 | 60 | 3-115-3067 | 0.977 | 0.960 | 0.921 | 0.901 |
| bsds388067_k8 | 80x53 | 22044 | 10.40 | 0.526 | 203 | 3-57-847 | 0.981 | 0.922 | 0.932 | 0.845 |
| bsds388067_k16 | 80x53 | 43262 | 20.41 | 0.580 | 203 | 3-57-847 | 0.967 | 0.870 | 0.910 | 0.783 |

## Caveats

- The ground truth is very coarse for a 4240-node graph: annotators draw
  2-92 segments, mostly 5-30, so groups are big (median group sizes 36-917 nodes,
  largest up to 89% of the image) and few. Any purely local graph puts most
  edges inside a group; the control already reaches 0.98 / 0.96 "within" at
  k = 8 / 16. The colour signal is the gap between graph and control, which
  is small: 0-6 points (k = 8) and 0.4-10 points (k = 16) on the pooled measure, and on one image (226022,
  k = 8) the colour graph is slightly *worse* than the grid on the
  annotator-mean measure.
- Pooling the annotators makes heavy overlap of near-copies: the same object
  outlined by different people gives groups that differ by a few blocks
  (13-115 pairs per image with Jaccard >= 0.9). Only exact duplicates were
  merged. A method that finds one copy is scored against all of them.
- Segments are not always spatially connected after the majority mapping,
  and segments thinner than about 3 blocks (18 pixels) vanish or shrink;
  1-35 segments per image fell under 3 nodes and were dropped.
- Colour similarity is a weak and sometimes wrong cue for object
  segmentation (texture, shading, objects with several colours, similar
  colours across a border).
- The graphs are planar-like local grids: n_components = 1 everywhere, and
  all edges span at most 3 blocks.
- `nodes.tsv` stores the mean colour of each block, i.e. an 80x53 thumbnail
  of a Corel photograph. The BSDS terms allow non-commercial research use;
  consider this before redistributing the folders outside the research
  context (the edges and groups alone do not need the colours).
- One dataset, one block size (6), one radius (3), one colour space, 15
  images.
