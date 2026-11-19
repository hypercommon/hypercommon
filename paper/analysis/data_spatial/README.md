# Visium spot graphs of human DLPFC with cortical layers as ground truth

60 graph folders from the 12 human dorsolateral prefrontal cortex (DLPFC) 10x
Visium slices of Maynard et al. 2021, with the manual layer annotation
(L1-L6, white matter WM) as ground truth. Built by
`paper/analysis/spatial_build.py` (run with `python -I`).

| set | folders | groups |
|---|---|---|
| `dlpfc_<S>_k8`, `dlpfc_<S>_k16` (12 slices) | 24 | one group per layer |
| `dlpfc_<S>_k8_null`, `dlpfc_<S>_k16_null` | 24 | the same |
| `dlpfc_<S>_parts_k8` (6 slices, see below) | 6 | each layer split into its spatially connected parts |
| `dlpfc_<S>_parts_k8_null` | 6 | the same |

`<S>` is the slice id (151507 ... 151676). Every `_null` folder is the
grid-only control for the folder of the same name without `_null`: same
nodes, same `groups.txt`, same `nodes.tsv`, different edges.

## Sources

All files are from the authors' own public locations; no account or login was
used and nothing downloaded was executed.

| what | URL |
|---|---|
| Space Ranger filtered count matrices (10x HDF5), one per slice | `https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/<S>_filtered_feature_bc_matrix.h5` |
| Spot positions (Space Ranger `tissue_positions_list.txt`), one per slice | `https://raw.githubusercontent.com/LieberInstitute/HumanPilot/044446d6bd8fc154aa74f7be62ec67effb1ec376/10X/<S>/tissue_positions_list.txt` |
| Spot-level layer annotation (`barcode`, slice, layer) | `https://raw.githubusercontent.com/LieberInstitute/HumanPilot/044446d6bd8fc154aa74f7be62ec67effb1ec376/10X/barcode_level_layer_map.tsv` |
| Per-slice Space Ranger metrics (used only for the checks below) | `https://raw.githubusercontent.com/LieberInstitute/HumanPilot/044446d6bd8fc154aa74f7be62ec67effb1ec376/Analysis/visium_dlpfc_pilot_sample_metrics.tsv` |
| Repository (official LIBD code for the paper) | https://github.com/LieberInstitute/HumanPilot |
| spatialLIBD (official R/Bioconductor data package) | https://bioconductor.org/packages/spatialLIBD, http://research.libd.org/spatialLIBD/ |

- Dataset: Maynard, Collado-Torres, Weber, Uytingco, Barry, Williams, Catallini,
  Tran, Besich, Tippani, Chew, Yin, Kleinman, Hyde, Rao, Hicks, Martinowich,
  Jaffe, "Transcriptome-scale spatial gene expression in the human dorsolateral
  prefrontal cortex", *Nature Neuroscience* 24:425-436 (2021),
  doi:10.1038/s41593-020-00787-0. Layer annotation: manual, by the authors,
  from cytoarchitecture and marker genes.
- The S3 links are the ones the HumanPilot repository lists in
  `AWS_File_locations.tsv` ("the raw data we received from 10x Genomics",
  README). The S3 objects are dated 2 March 2020. GitHub files are pinned to
  commit `044446d6bd8f` (the newest on 13 November 2026).
- R was not used. The spatialLIBD `.rds`/ExperimentHub objects were not read,
  and the third-party Hugging Face repackaging was not used.
- Checks against the official numbers: for every slice the filtered matrix has
  exactly as many spots as Space Ranger's "Estimated Number of Spots" in the
  metrics table, and all of them are `in_tissue = 1` in the positions file;
  the 12 slices hold 47,681 spots, the number in the paper and in the
  spatialLIBD object. 47,329 spots carry a layer label; 352 have none and
  are dropped. Slice 151673 has 3,611 labelled spots (L1 273, L2 253, L3 989,
  L4 218, L5 673, L6 692, WM 513), the counts used by the usual DLPFC
  benchmarks. Every labelled barcode is in its slice's filtered matrix.
- Download SHA-256 (first 12 hex digits): layer map `b5d81f0e6181`; h5
  151507 `686521ab02a6`, 151508 `48ab74b340c7`, 151509 `35a4850aebdf`, 151510
  `421a77bad2e0`, 151669 `4c92c56291f2`, 151670 `931f50b5de5a`, 151671
  `6b6b86696b37`, 151672 `a28ad2b88314`, 151673 `216e8010e1c3`, 151674
  `0ecbabdb0ea9`, 151675 `eba1f479b85d`, 151676 `f436605010a0`.

Licence / terms:
- The HumanPilot repository and the spatialLIBD package are released under
  the Artistic License 2.0 (repository `LICENSE` file; spatialLIBD
  `DESCRIPTION`). No separate licence is stated for the data files; they are
  posted for public download without login by the authors. Cite Maynard et
  al. 2021 (and spatialLIBD, Pardo et al. 2022, *BMC Genomics* 23:434) when
  using these graphs.

## Construction (rule fixed before any detection was run)

1. **Nodes.** The spots of the slice's filtered matrix that carry a layer
   label. Ordered by `(array_row, array_col)` and numbered 0..n-1;
   `nodes.tsv` gives barcode, grid position, full-resolution pixel position and
   layer.
2. **Expression features** (on the n nodes only). Genes with a count > 0 in
   at least 3 nodes; counts scaled to 10,000 per spot, then `log1p`; the 2,000
   most variable genes by the Seurat v1 dispersion rule (scanpy
   `flavor="seurat"`: dispersion of the un-logged normalised data, z-scored
   within 20 equal-width bins of log mean); each gene centred and scaled to
   unit variance, clipped at 10; exact PCA to 30 components. No spatial
   smoothing, no batch correction, every slice processed alone.
3. **Candidates.** For node i, the nodes at hexagonal grid distance 1, 2 or 3
   (the first three rings, at most 6 + 12 + 18 = 36 spots; fewer at the tissue
   edge). On the Visium grid the six neighbours of `(r, c)` are `(r, c+-2)`
   and `(r+-1, c+-1)`; grid distance = `|dr| + max(0, (|dc| - |dr|) / 2)`.
4. **Edges (exact rule).** Node i keeps the k candidates with the smallest
   Euclidean distance to it in the 30-dimensional PCA space (ties: smaller
   centre-to-centre distance on the grid, then smaller id). The undirected,
   unweighted edge set is the union of these choices (so degree >= k, except
   at the tissue edge where a spot can have fewer than k candidates). k = 8
   and 16. No self-loops, no duplicates, each edge written once as
   `u<TAB>v` with u < v.
5. **Control (`_null`).** The same rule with the expression score replaced
   by the grid: node i keeps its k candidates with the smallest
   centre-to-centre distance (spot pitch units: 6 at 1, 6 at sqrt(3), 6 at 2,
   then 12 at sqrt(7), ...), ties to the smaller id; union. So the control
   is a pure triangular-lattice graph with no expression signal, on exactly
   the same nodes and candidate sets.
6. **Groups.**
   - *layers* (all slices, both k): one group per layer present in the slice
     (7 on 8 slices; 5, L3-WM, on 151669-151672, which have no L1/L2).
   - *parts* (k = 8 only): each layer split into the connected components of
     the candidate graph (grid distance <= 3) restricted to that layer's
     nodes; parts with fewer than 3 nodes are dropped (none were). On 6 slices
     the parts are exactly the layers, so a parts set was written only for the
     6 where they differ: 151507, 151508, 151509, 151510 (L2 and L3 each in
     two parts), 151673 (L1 in three parts: 232, 38, 3 spots), 151675 (L1 in
     two: 290, 38).

Why both group sets, and why layers first. Layers are the standard ground
truth for this benchmark (every spatial-domain paper scores against them), so
they are the main set. But a layer that appears as two separate bands on the
slide cannot be one community of any local graph method: no edge (and no
candidate pair) joins the two bands. The parts set is the ground truth a
purely local method could reach. Because every within-layer edge is
necessarily within a part, the edge counts and `frac_edges_within_group` are
the same in both sets; only the group lists differ. Parts were made only at
k = 8 to stay within the folder budget (60). Both sets are partitions: every
node is in exactly one group, and groups are large (166-2,175 nodes), not the
small groups of the Hi-C and protein families.

Why these choices: the candidate radius keeps every edge local (at most 3
spot pitches, about 300 um) so that, as in the Hi-C top-k recipe, the graph is
a lattice whose links are chosen by the data; k = 8 and 16 bracket the 6 and
18 spots in the first one and two rings. Euclidean distance on 30 PCs of
log-normalised HVG expression is the default of the common pipelines
(Seurat/scanpy) and was taken without trying alternatives.

## Reproduce

Each download in its own folder under one RAW directory:
`RAW/h5_<S>/<S>_filtered_feature_bc_matrix.h5`,
`RAW/pos_<S>/tissue_positions_list.txt`,
`RAW/layer_map/barcode_level_layer_map.tsv`. Then, from `paper/`:

    python -I analysis/spatial_build.py inspect --raw RAW      # parts and candidate counts only
    python -I analysis/spatial_build.py build   --raw RAW --out analysis/data_spatial

`build` writes all 60 folders and `summary.json` (about 11 min; the defaults
are `--ks 8 16 --groups both --parts-ks 8`). The controls are written by the
same command. Built with Python 3.12.10, numpy 2.5.4, scipy 1.18.1, h5py
3.16.0, networkx 3.7, scikit-learn 1.9.1 (one diagnostic only). The script
sets `OPENBLAS_NUM_THREADS=1` (multi-threaded BLAS was ~30x slower on this
machine); results are deterministic.

## Files

Each folder has:

- `edges.tsv`: one undirected edge per line, `u<TAB>v`, u < v, ids 0..n-1, no header.
- `groups.txt`: one group per line, space-separated sorted ids; layers in the
  order L1 ... L6, WM (parts: by layer, then larger part first).
- `nodes.tsv`: header `id barcode array_row array_col pxl_row_in_fullres pxl_col_in_fullres layer`.

`summary.json` has one row per folder (`control` true/false). Rows for the
expression graphs also carry the numbers of their control (`null_*`) and
diagnostics: share of edges at grid distance 1/2/3, share of edges shared
with the control, the same top-k rule with a seeded random score
(`random_rule_*`), the share of all candidate pairs inside a group
(`candidate_pairs_frac_within_group`), and an expression kNN graph with no
spatial restriction (`global_expr_knn_*`). The diagnostic graphs were not
written as folders.

## Counts

Per slice. "dropped" = in-tissue spots without a layer label.

| slice | spots in tissue | labelled (n) | dropped | layers | layer sizes (L1 L2 L3 L4 L5 L6 WM) | candidates min / median | parts differ |
|---|---|---|---|---|---|---|---|
| 151507 | 4226 | 4221 | 5 | 7 | 817 305 1215 369 675 486 354 | 11 / 36 | yes |
| 151508 | 4384 | 4381 | 3 | 7 | 866 295 1385 373 737 525 200 | 8 / 36 | yes |
| 151509 | 4789 | 4788 | 1 | 7 | 1189 602 1884 369 363 215 166 | 11 / 36 | yes |
| 151510 | 4634 | 4595 | 39 | 7 | 1180 650 1774 318 310 179 184 | 11 / 36 | yes |
| 151669 | 3661 | 3636 | 25 | 5 | - - 2141 364 510 391 230 | 7 / 36 | no |
| 151670 | 3498 | 3484 | 14 | 5 | - - 2175 211 581 308 209 | 8 / 36 | no |
| 151671 | 4110 | 4093 | 17 | 5 | - - 1918 245 721 760 449 | 7 / 36 | no |
| 151672 | 4015 | 3888 | 127 | 5 | - - 1575 304 728 882 399 | 10 / 36 | no |
| 151673 | 3639 | 3611 | 28 | 7 | 273 253 989 218 673 692 513 | 1 / 36 | yes |
| 151674 | 3673 | 3635 | 38 | 7 | 380 224 924 247 621 614 625 | 10 / 36 | no |
| 151675 | 3592 | 3566 | 26 | 7 | 328 275 771 275 732 533 652 | 9 / 36 | yes |
| 151676 | 3460 | 3431 | 29 | 7 | 289 254 836 254 649 616 533 | 6 / 36 | no |

Per graph. "within" = fraction of edges whose two ends are in the same group.
"random rule" = the same top-k rule with a random score (diagnostic). All
graphs and controls are connected except slice 151673 (2 components in both
the graph and its control: a 3-spot island more than 3 pitches from every
other labelled spot).

| graph | n | m | mean degree | clustering | components | groups | group size min-median-max | within | within, control | within, random rule | m control | clustering control | edges at grid distance 1 / 2 / 3 | edges shared with control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dlpfc_151507_k8 | 4221 | 25704 | 12.18 | 0.325 | 1 | 7 | 305-486-1215 | 0.846 | 0.914 | 0.841 | 21053 | 0.419 | 0.20 / 0.34 / 0.46 | 0.32 |
| dlpfc_151507_parts_k8 | 4221 | 25704 | 12.18 | 0.325 | 1 | 9 | 73-369-1142 | 0.846 | 0.914 | 0.841 | 21053 | 0.419 | 0.20 / 0.34 / 0.46 | 0.32 |
| dlpfc_151507_k16 | 4221 | 46441 | 22.00 | 0.435 | 1 | 7 | 305-486-1215 | 0.845 | 0.882 | 0.841 | 34517 | 0.520 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151508_k8 | 4381 | 26993 | 12.32 | 0.324 | 1 | 7 | 200-525-1385 | 0.851 | 0.919 | 0.849 | 21866 | 0.421 | 0.19 / 0.35 / 0.47 | 0.32 |
| dlpfc_151508_parts_k8 | 4381 | 26993 | 12.32 | 0.324 | 1 | 9 | 88-373-1259 | 0.851 | 0.919 | 0.849 | 21866 | 0.421 | 0.19 / 0.35 / 0.47 | 0.32 |
| dlpfc_151508_k16 | 4381 | 48620 | 22.20 | 0.437 | 1 | 7 | 200-525-1385 | 0.853 | 0.888 | 0.849 | 35824 | 0.520 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151509_k8 | 4788 | 29014 | 12.12 | 0.323 | 1 | 7 | 166-369-1884 | 0.861 | 0.917 | 0.848 | 23887 | 0.419 | 0.20 / 0.34 / 0.46 | 0.33 |
| dlpfc_151509_parts_k8 | 4788 | 29014 | 12.12 | 0.323 | 1 | 9 | 166-363-1231 | 0.861 | 0.917 | 0.848 | 23887 | 0.419 | 0.20 / 0.34 / 0.46 | 0.33 |
| dlpfc_151509_k16 | 4788 | 52462 | 21.91 | 0.435 | 1 | 7 | 166-369-1884 | 0.855 | 0.887 | 0.847 | 39056 | 0.518 | 0.19 / 0.34 / 0.47 | 0.50 |
| dlpfc_151510_k8 | 4595 | 27910 | 12.15 | 0.322 | 1 | 7 | 179-318-1774 | 0.848 | 0.916 | 0.844 | 22944 | 0.420 | 0.19 / 0.34 / 0.46 | 0.32 |
| dlpfc_151510_parts_k8 | 4595 | 27910 | 12.15 | 0.322 | 1 | 9 | 179-318-1180 | 0.848 | 0.916 | 0.844 | 22944 | 0.420 | 0.19 / 0.34 / 0.46 | 0.32 |
| dlpfc_151510_k16 | 4595 | 50515 | 21.99 | 0.432 | 1 | 7 | 179-318-1774 | 0.846 | 0.884 | 0.844 | 37522 | 0.519 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151669_k8 | 3636 | 22145 | 12.18 | 0.322 | 1 | 5 | 230-391-2141 | 0.909 | 0.941 | 0.903 | 18113 | 0.426 | 0.19 / 0.33 / 0.47 | 0.32 |
| dlpfc_151669_k16 | 3636 | 39787 | 21.89 | 0.436 | 1 | 5 | 230-391-2141 | 0.906 | 0.926 | 0.901 | 29817 | 0.525 | 0.18 / 0.34 / 0.48 | 0.51 |
| dlpfc_151670_k8 | 3484 | 21164 | 12.15 | 0.316 | 1 | 5 | 209-308-2175 | 0.906 | 0.940 | 0.902 | 17357 | 0.424 | 0.19 / 0.34 / 0.47 | 0.32 |
| dlpfc_151670_k16 | 3484 | 38184 | 21.92 | 0.434 | 1 | 5 | 209-308-2175 | 0.904 | 0.925 | 0.902 | 28523 | 0.524 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151671_k8 | 4093 | 25048 | 12.24 | 0.319 | 1 | 5 | 245-721-1918 | 0.902 | 0.941 | 0.900 | 20412 | 0.422 | 0.19 / 0.34 / 0.47 | 0.32 |
| dlpfc_151671_k16 | 4093 | 45180 | 22.08 | 0.435 | 1 | 5 | 245-721-1918 | 0.901 | 0.926 | 0.899 | 33452 | 0.521 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151672_k8 | 3888 | 23896 | 12.29 | 0.321 | 1 | 5 | 304-728-1575 | 0.899 | 0.940 | 0.896 | 19382 | 0.422 | 0.19 / 0.34 / 0.47 | 0.32 |
| dlpfc_151672_k16 | 3888 | 42970 | 22.10 | 0.434 | 1 | 5 | 304-728-1575 | 0.897 | 0.924 | 0.896 | 31826 | 0.522 | 0.18 / 0.34 / 0.48 | 0.50 |
| dlpfc_151673_k8 | 3611 | 21774 | 12.06 | 0.329 | 2 | 7 | 218-513-989 | 0.813 | 0.875 | 0.805 | 17978 | 0.424 | 0.20 / 0.35 / 0.45 | 0.33 |
| dlpfc_151673_parts_k8 | 3611 | 21774 | 12.06 | 0.329 | 2 | 9 | 3-253-989 | 0.813 | 0.875 | 0.805 | 17978 | 0.424 | 0.20 / 0.35 / 0.45 | 0.33 |
| dlpfc_151673_k16 | 3611 | 39410 | 21.83 | 0.444 | 2 | 7 | 218-513-989 | 0.806 | 0.854 | 0.803 | 29598 | 0.524 | 0.19 / 0.34 / 0.47 | 0.51 |
| dlpfc_151674_k8 | 3635 | 21696 | 11.94 | 0.348 | 1 | 7 | 224-614-924 | 0.818 | 0.874 | 0.802 | 18095 | 0.425 | 0.22 / 0.35 / 0.43 | 0.35 |
| dlpfc_151674_k16 | 3635 | 39413 | 21.69 | 0.456 | 1 | 7 | 224-614-924 | 0.809 | 0.852 | 0.800 | 29895 | 0.525 | 0.20 / 0.34 / 0.46 | 0.53 |
| dlpfc_151675_k8 | 3566 | 21782 | 12.22 | 0.335 | 1 | 7 | 275-533-771 | 0.807 | 0.871 | 0.793 | 17779 | 0.427 | 0.20 / 0.34 / 0.46 | 0.33 |
| dlpfc_151675_parts_k8 | 3566 | 21782 | 12.22 | 0.335 | 1 | 8 | 38-412-771 | 0.807 | 0.871 | 0.793 | 17779 | 0.427 | 0.20 / 0.34 / 0.46 | 0.33 |
| dlpfc_151675_k16 | 3566 | 39170 | 21.97 | 0.447 | 1 | 7 | 275-533-771 | 0.801 | 0.851 | 0.791 | 29332 | 0.527 | 0.19 / 0.34 / 0.47 | 0.52 |
| dlpfc_151676_k8 | 3431 | 20819 | 12.14 | 0.328 | 1 | 7 | 254-533-836 | 0.801 | 0.872 | 0.794 | 17091 | 0.427 | 0.19 / 0.34 / 0.47 | 0.32 |
| dlpfc_151676_k16 | 3431 | 37589 | 21.91 | 0.443 | 1 | 7 | 254-533-836 | 0.798 | 0.849 | 0.794 | 28164 | 0.526 | 0.18 / 0.34 / 0.48 | 0.51 |

Medians over the 12 slices (layers): within 0.850 (k = 8) and 0.850
(k = 16); control 0.917 and 0.886; random rule 0.846 and 0.846; all
candidate pairs 0.845. Clustering 0.32 / 0.44 (graph) against 0.42 / 0.52
(control). Clustering is networkx `average_clustering`.

## Caveats

- **The expression signal at spot level is weak; the graphs are close to
  random geometric graphs of radius 3.** Within the 36 candidates, the
  expression rule puts an edge inside a layer barely more often than a random
  choice (median 0.850 vs 0.846) and clearly less often than the grid control
  (0.917 at k = 8, 0.886 at k = 16), on every slice. Its edges spread over
  the three rings roughly in proportion to their size (about 0.19 / 0.34 /
  0.47, against 6 / 12 / 18 of 36 candidates), only a third (k = 8) to a half
  (k = 16) of them are control edges, and clustering is lower than the
  control's. A kNN graph on the same PCA space without any spatial limit has
  only 40-62% of its edges inside a layer (`global_expr_knn_*`). This fits the
  known picture that single Visium spots are noisy and layers are told apart
  mainly by spatial context (methods that smooth over neighbours do much
  better than expression-only clustering). The rule was fixed before these
  numbers were seen and was not changed after; a smaller radius or spatially
  smoothed features would move the graph towards the control, but would also
  be a choice made with knowledge of the ground truth.
- So the high "within" values (0.80-0.91) come almost entirely from the
  geometry: layers are bands tens of spots wide, and 79-90% of all candidate
  pairs are already inside one layer. For HC-style detection the groups are
  few (5-9) and very large, so the task is mostly to find a handful of long
  borders in a 2D sheet, not many small dense groups.
- The control is a genuine lattice (triangular grid, k nearest), so on these
  slices the control is more "lattice-like" than the data graph, the
  reverse of the intended contrast.
- Layer boundaries are the authors' manual annotation;
  borders between neighbouring layers are uncertain by a few spots, and
  L2/L3 and L5/L6 differ only gradually in expression.
- 151669-151672 (one donor) have no L1/L2 annotation. The 12 slices come
  from 3 donors, 4 sections each (two pairs of adjacent 10 um sections, the
  pairs 300 um apart), so the 12 graphs are not independent.
- At the tissue edge and next to unlabelled spots a node has fewer than 36
  candidates (median 36; fewer than k = 16 candidates for 7-35 nodes per
  slice), so a few nodes there have degree < k.
