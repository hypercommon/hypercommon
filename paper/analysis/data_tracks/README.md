# Detection graphs from pedestrian tracking with track IDs as ground truth

30 graphs: 15 frame windows from the MOT17 training sequences, each with
k = 8 and k = 16, plus a distance-only control `<name>_null` for every graph
(60 folders in all). Nodes are ground-truth pedestrian boxes, edges join a box
to the boxes nearby in space and time that overlap it most, and the groups are
the track IDs. Each group is a long, thin chain through time (one box per
frame). Built by `paper/analysis/tracks_build.py` (run with `python -I`).

## Sources

| What | URL |
|---|---|
| MOT17 labels archive (gt + public detections, no images) | https://motchallenge.net/data/MOT17Labels.zip |
| Dataset page (licence, sequence table) | https://motchallenge.net/data/MOT17/ |

- Dataset: A. Milan, L. Leal-Taixe, I. Reid, S. Roth, K. Schindler,
  "MOT16: A Benchmark for Multi-Object Tracking", arXiv:1603.00831 (2016).
  MOT17 uses the MOT16 sequences with a new, more accurate ground truth.
- Downloaded on 2026-11-13 with curl, no login, no account.
- Archive: 10,107,022 bytes,
  sha256 `0aa79322e91583369f42f17c4d79a0b145380d8732487bba59272048dc82b2b9`,
  md5 `26048d544834c59e73b80b796c392ab3`. MOTChallenge publishes no
  checksum. Only `train/*/gt/gt.txt` and `train/*/seqinfo.ini` were
  extracted and read (numpy `loadtxt`, `configparser`). No file from the
  archive was executed. Each sequence comes in three copies (DPM, FRCNN, SDP
  detections); their `gt.txt` files are byte-identical (checked by md5), so
  the FRCNN copy is used. The detector outputs (`det/`) are not used.
- Licence: the MOT17 page says the data are published under
  **CC BY-NC-SA 3.0** (attribution, non-commercial, share-alike), and that
  for sequences taken from third-party sources the original conditions may
  apply too (here MOT17-05, from Ess et al. ICCV 2007, and MOT17-09/-11 from
  MOT16). These graphs are derived from the labels, so they fall under the
  same licence: cite Milan et al. 2016, no commercial use, share alike.

## Selection (fixed before building)

Per sequence, N = number of kept boxes (see Construction 1) and
c = floor(N / 3500 + 0.5) windows. Frames 1..seqLength are cut into c runs of
consecutive frames; cut j (j = 1..c-1) goes after the first frame where the
running box count reaches j * N / c. This gave 24 windows:

| sequence | fps | frames | rows in gt.txt | kept boxes | windows |
|---|---|---|---|---|---|
| MOT17-02 | 30 | 600 | 30,003 | 9,393 | 3 |
| MOT17-04 | 30 | 1050 | 108,005 | 39,821 | 11 |
| MOT17-05 | 14 | 837 | 8,013 | 4,767 | 1 |
| MOT17-09 | 30 | 525 | 10,411 | 3,745 | 1 |
| MOT17-10 | 30 | 654 | 17,450 | 10,560 | 3 |
| MOT17-11 | 30 | 900 | 10,617 | 7,518 | 2 |
| MOT17-13 | 25 | 750 | 20,202 | 10,024 | 3 |

With `rng = random.Random(2026)`: first one window per sequence,
`rng.choice` over that sequence's windows, sequences in sorted order; then
`rng.sample` of the 17 remaining windows (sorted) to reach 15 in all. Chosen
(window number is 1-based inside its sequence):

```
02 w1 w2 w3 | 04 w2 w3 w6 w9 w11 | 05 w1 | 09 w1 | 10 w3 | 11 w1 w2 | 13 w1 w3
```

Every training sequence is used. No window was swapped or re-drawn. 15
windows (not 20) because two k values plus a control per graph would
otherwise exceed the ~60-folder budget of the downstream benchmark.

## Construction

All rules were fixed before any graph was built. No community detection was
run, and nothing was tuned; only the structural statistics below were
computed.

1. **Nodes.** Rows of `gt.txt` (`frame, id, left, top, w, h, flag, class,
   visibility`) with flag = 1 (counted in evaluation), class = 1
   (pedestrian) and visibility >= 0.25. One node per box. Inside a window the
   boxes are sorted by (frame, left, top, w, h, line in gt.txt) and numbered
   0..n-1; the track id plays no part in the order. Each track has at most
   one box per frame (checked).
2. **Candidates.** For node a: every other node b with
   |frame_b - frame_a| <= 3 (same frame included) and centre distance
   <= 1.5 x h_a (centre = (left + w/2, top + h/2), h_a = a's box height).
   Not symmetric (it uses a's height).
3. **Edges (exact rule).** Node a keeps the k candidates with the largest
   IoU of the two boxes as given (no motion shift); ties by smaller centre
   distance, then smaller id; all candidates if there are fewer than k. The
   undirected, unweighted edge set is the union over a. No self-loops, each
   edge once as `u<TAB>v` with u < v.
4. **Groups.** One group per track id with >= 3 nodes in the window. Groups
   are disjoint. Nodes of shorter tracks (0-8 tracks per window) stay in the
   graph but in no group; 99.8-100% of nodes are in a group.
5. **Control `<name>_null`.** Same nodes, same candidate sets, same
   `groups.txt` and `nodes.tsv`. Node a keeps the k candidates with the
   smallest centre distance (ties by smaller id); union over a. No overlap
   signal, only distance.

Why this rule: it is the Hi-C / image top-k recipe on a space-time order. The
radius keeps every edge local; choosing by box overlap inside the radius
should favour the same person in the next frames. The +-3 frame window and
1.5 x height radius were set from what a person can move in 0.1 s, not from
any graph statistic.

## Reproduce

Needs numpy and scipy (we used a separate venv: Python 3.12, numpy 2.5.4,
scipy 1.18.1). Put the archive in its own empty folder and extract only the
labels there (`unzip MOT17Labels.zip 'train/*/gt/gt.txt' 'train/*/seqinfo.ini' -d ex`),
then from `paper/`:

    python -I analysis/tracks_build.py build --mot <folder holding train/> --out analysis/data_tracks

It writes all 60 folders and `summary.json` in under a minute on one core.
Two runs give byte-identical output.

## Files

Each folder `mot17_<SEQ>_w<J>_k<K>/` and `mot17_<SEQ>_w<J>_k<K>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v.
- `groups.txt` - one ground-truth group (one track) per line,
  space-separated sorted node ids, size >= 3.
- `nodes.tsv` - header `id frame x y w h`: the original frame number in the
  sequence and the box (`x, y` = top-left corner in pixels, as in gt.txt).
  The track id is not stored.

`summary.json` has one row per graph, with the control's numbers under
`null_*`. Fields: sequence, window, frame range, fps, image size, k, n, m,
mean/min/max degree, clustering (average local clustering, 0 for degree < 2,
same as networkx), n_components, candidates per node (mean, min) and the
fraction of nodes with fewer than k candidates, n_tracks, n_groups,
n_tracks_dropped_small, group size min/mean/median/max, group span in frames
(last - first + 1; mean, max), frac_nodes_in_group, and:

- `frac_edges_within_group`: fraction of edges whose two ends are in the same
  group (the "within" measure used for the other families);
- `n_within_candidate_pairs`: unordered same-group pairs where one end is a
  candidate of the other;
- `within_ceiling`: n_within_candidate_pairs / m (the largest "within" an
  edge set of this size could have under the candidate rule);
- `frac_within_candidate_pairs_linked`: share of those pairs that became an
  edge;
- `frac_groups_connected`, `group_pieces_mean`: whether each track's
  subgraph is connected, and its mean number of connected pieces;
- `frac_edges_same_frame`, `mean_edge_dt`: time span of the edges.

## Counts

n = 3123-4767 (53,684 boxes in the 15 windows). Groups per window 26-115,
mean track length 41-144 nodes (= frames with a kept box), mean span 50-177
frames. Median over the 15 windows (min-max in brackets):

| quantity | k = 8 | k = 8 null | k = 16 | k = 16 null |
|---|---|---|---|---|
| m | 16,593 (14,190-22,006) | 16,529 (14,132-22,236) | 33,222 (26,581-44,009) | 32,728 (26,493-45,882) |
| mean degree | 9.23 | 9.22 | 18.49 | 18.37 |
| clustering | 0.453 (0.433-0.517) | 0.451 (0.441-0.524) | 0.577 (0.558-0.605) | 0.578 (0.556-0.605) |
| within | 0.626 (0.564-0.670) | 0.630 (0.577-0.664) | 0.311 (0.296-0.349) | 0.318 (0.284-0.354) |
| within ceiling | 0.629 (0.592-0.673) | | 0.311 (0.296-0.359) | |
| same-group candidate pairs linked | 1.000 (0.838-1.000) | 1.000 (0.861-1.000) | 1.000 (0.971-1.000) | 1.000 (0.982-1.000) |
| pieces per track | 1.50 (1.05-2.38) | same | 1.50 (1.05-2.38) | same |

Per graph (within = `frac_edges_within_group`; linked =
`frac_within_candidate_pairs_linked`):

| graph | frames | n | m | mean degree | clustering | groups | size min-mean-max | span mean | pieces/track | within | within, null | ceiling | linked | linked, null |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mot17_02_w1_k8 | 1-207 | 3133 | 14536 | 9.28 | 0.453 | 28 | 4-111.9-207 | 131.6 | 1.50 | 0.629 | 0.622 | 0.629 | 1.000 | 1.000 |
| mot17_02_w1_k16 | 1-207 | 3133 | 28355 | 18.10 | 0.574 | 28 | 4-111.9-207 | 131.6 | 1.50 | 0.322 | 0.321 | 0.322 | 1.000 | 1.000 |
| mot17_02_w2_k8 | 208-398 | 3137 | 14399 | 9.18 | 0.478 | 38 | 5-82.6-191 | 104.2 | 1.76 | 0.626 | 0.626 | 0.626 | 1.000 | 1.000 |
| mot17_02_w2_k16 | 208-398 | 3137 | 29007 | 18.49 | 0.568 | 38 | 5-82.6-191 | 104.2 | 1.76 | 0.311 | 0.309 | 0.311 | 1.000 | 1.000 |
| mot17_02_w3_k8 | 399-600 | 3123 | 14541 | 9.31 | 0.464 | 37 | 3-84.4-194 | 120.9 | 2.14 | 0.611 | 0.607 | 0.611 | 1.000 | 1.000 |
| mot17_02_w3_k16 | 399-600 | 3123 | 29067 | 18.61 | 0.573 | 37 | 3-84.4-194 | 120.9 | 2.14 | 0.306 | 0.302 | 0.306 | 1.000 | 1.000 |
| mot17_04_w2_k8 | 107-206 | 3593 | 16593 | 9.24 | 0.448 | 39 | 7-92.1-100 | 92.7 | 1.08 | 0.635 | 0.637 | 0.635 | 1.000 | 1.000 |
| mot17_04_w2_k16 | 107-206 | 3593 | 32535 | 18.11 | 0.558 | 39 | 7-92.1-100 | 92.7 | 1.08 | 0.324 | 0.323 | 0.324 | 1.000 | 1.000 |
| mot17_04_w3_k8 | 207-299 | 3620 | 16905 | 9.34 | 0.440 | 43 | 11-84.2-93 | 85.8 | 1.12 | 0.625 | 0.630 | 0.625 | 1.000 | 1.000 |
| mot17_04_w3_k16 | 207-299 | 3620 | 33507 | 18.51 | 0.563 | 43 | 11-84.2-93 | 85.8 | 1.12 | 0.316 | 0.314 | 0.316 | 1.000 | 1.000 |
| mot17_04_w6_k8 | 479-571 | 3632 | 17075 | 9.40 | 0.433 | 43 | 20-84.5-93 | 85.6 | 1.05 | 0.622 | 0.630 | 0.622 | 1.000 | 1.000 |
| mot17_04_w6_k16 | 479-571 | 3632 | 35306 | 19.44 | 0.581 | 43 | 20-84.5-93 | 85.6 | 1.05 | 0.301 | 0.303 | 0.301 | 1.000 | 1.000 |
| mot17_04_w9_k8 | 776-874 | 3613 | 17140 | 9.49 | 0.445 | 43 | 11-84.0-99 | 86.5 | 1.09 | 0.616 | 0.616 | 0.616 | 1.000 | 1.000 |
| mot17_04_w9_k16 | 776-874 | 3613 | 34407 | 19.05 | 0.558 | 43 | 11-84.0-99 | 86.5 | 1.09 | 0.307 | 0.318 | 0.307 | 1.000 | 1.000 |
| mot17_04_w11_k8 | 967-1050 | 3587 | 16980 | 9.47 | 0.438 | 49 | 5-73.2-84 | 74.2 | 1.06 | 0.615 | 0.632 | 0.615 | 1.000 | 1.000 |
| mot17_04_w11_k16 | 967-1050 | 3587 | 33846 | 18.87 | 0.562 | 49 | 5-73.2-84 | 74.2 | 1.06 | 0.308 | 0.319 | 0.308 | 1.000 | 1.000 |
| mot17_05_w1_k8 | 1-837 | 4767 | 22006 | 9.23 | 0.496 | 115 | 3-41.4-376 | 50.4 | 1.71 | 0.584 | 0.577 | 0.592 | 0.986 | 0.985 |
| mot17_05_w1_k16 | 1-837 | 4767 | 44009 | 18.46 | 0.598 | 115 | 3-41.4-376 | 50.4 | 1.70 | 0.296 | 0.284 | 0.296 | 1.000 | 1.000 |
| mot17_09_w1_k8 | 1-525 | 3745 | 16784 | 8.96 | 0.450 | 26 | 5-144.0-327 | 177.1 | 2.38 | 0.647 | 0.646 | 0.647 | 1.000 | 1.000 |
| mot17_09_w1_k16 | 1-525 | 3745 | 34048 | 18.18 | 0.592 | 26 | 5-144.0-327 | 177.1 | 2.38 | 0.319 | 0.318 | 0.319 | 1.000 | 1.000 |
| mot17_10_w3_k8 | 429-654 | 3515 | 15883 | 9.04 | 0.458 | 28 | 18-125.5-226 | 134.2 | 1.50 | 0.646 | 0.640 | 0.648 | 0.998 | 0.997 |
| mot17_10_w3_k16 | 429-654 | 3515 | 33222 | 18.90 | 0.577 | 28 | 18-125.5-226 | 134.2 | 1.50 | 0.310 | 0.306 | 0.310 | 1.000 | 1.000 |
| mot17_11_w1_k8 | 1-408 | 3762 | 16309 | 8.67 | 0.471 | 34 | 3-110.5-408 | 119.7 | 1.41 | 0.670 | 0.664 | 0.673 | 0.996 | 0.995 |
| mot17_11_w1_k16 | 1-408 | 3762 | 31910 | 16.96 | 0.589 | 34 | 3-110.5-408 | 119.7 | 1.41 | 0.344 | 0.342 | 0.344 | 1.000 | 1.000 |
| mot17_11_w2_k8 | 409-900 | 3756 | 17083 | 9.10 | 0.453 | 44 | 4-85.3-459 | 99.4 | 1.50 | 0.635 | 0.623 | 0.636 | 0.998 | 0.998 |
| mot17_11_w2_k16 | 409-900 | 3756 | 35639 | 18.98 | 0.580 | 44 | 4-85.3-459 | 99.4 | 1.50 | 0.305 | 0.296 | 0.305 | 1.000 | 1.000 |
| mot17_13_w1_k8 | 1-143 | 3360 | 14642 | 8.72 | 0.472 | 46 | 7-73.0-143 | 80.1 | 1.63 | 0.633 | 0.636 | 0.656 | 0.966 | 0.974 |
| mot17_13_w1_k16 | 1-143 | 3360 | 28979 | 17.25 | 0.605 | 46 | 7-73.0-143 | 80.1 | 1.61 | 0.331 | 0.330 | 0.331 | 1.000 | 1.000 |
| mot17_13_w3_k8 | 350-750 | 3341 | 14190 | 8.49 | 0.517 | 49 | 5-68.1-238 | 73.6 | 1.67 | 0.564 | 0.582 | 0.673 | 0.838 | 0.861 |
| mot17_13_w3_k16 | 350-750 | 3341 | 26581 | 15.91 | 0.596 | 49 | 5-68.1-238 | 73.6 | 1.43 | 0.349 | 0.354 | 0.359 | 0.971 | 0.982 |

## Caveats

- **The overlap rule adds no group signal over distance.** Within +-3 frames
  a person's own boxes are almost always both the nearest and the most
  overlapping candidates, so both the IoU graph and the distance-only control
  link essentially every same-track candidate pair (median 100%). Their
  "within" values match to about 1 point (real minus control from -1.8 to
  +1.2 points; the real graph is ahead on 9 of 15 windows at k = 8 and on 11
  of 15 at k = 16). The two graphs differ only in which cross-track edges fill the
  remaining slots. The control is therefore not a weaker baseline here; it is
  about as good.
- **"Within" is capped by construction.** A node has at most 6 same-track
  candidates (one box per frame, +-3 frames), so at most 6 of its k choices
  can stay inside its group. "Within" sits at its ceiling: about 0.63 at
  k = 8 and 0.31 at k = 16. Most edges at k = 16 cross tracks
  (to people next to or crossing the box). This is far lower than the
  image family (0.97-0.99).
- **Chains are broken.** The visibility >= 0.25 filter removes occluded
  boxes, and a gap of more than 3 frames cuts a track's chain: tracks fall
  into 1.05-2.38 connected pieces on average (MOT17-04, few gaps: ~1.1;
  MOT17-09 and -02 w3, many occlusions: > 2). Only 35-95% of tracks are
  connected inside their own group. A method can only find a broken track
  as one group through the cross-track edges.
- Groups are long, thin chains: track lengths 3-459 nodes, and inside a
  track edges only span +-3 frames (a band). Same-frame edges (always
  cross-track) are 2-5% of edges at k = 8 and 8-10% at k = 16.
- Some nodes have no candidate at all (isolated people), so degree_min = 0
  and some graphs have up to 10 components (MOT17-13 w3).
- MOT17-05 (640x480, 14 fps, moving camera) has only one window (the whole
  sequence) and many short tracks (115 groups, mean 41 nodes); MOT17-04
  (elevated static camera, dense) gives 5 of the 15 windows.
- Windows are cut at fixed frames, so tracks running across a cut are
  truncated (their pieces in other windows are separate graphs).
- One dataset, one time radius (3), one space radius (1.5 x height), one
  visibility cut, 15 windows.
