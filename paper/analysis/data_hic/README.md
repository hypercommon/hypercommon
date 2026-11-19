# Hi-C chromosome graphs with Arrowhead contact domains as ground truth

18 graphs: GM12878 chromosomes 17, 18, 19, 20, 21, 22 at 25 kb, each with
k = 8, 16, 32. Built by `paper/analysis/hic_build.py`.

## Sources

| What | URL |
|---|---|
| Contact map (Juicer `.hic`, v7, hg19) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_insitu_primary%2Breplicate_combined_30.hic |
| Contact domains (Arrowhead) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary%2Breplicate_Arrowhead_domainlist.txt.gz |
| GEO series page | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525 |

- Dataset: Rao, Huntley, Durand, Stamenova, Bochkov, Robinson, Sanborn, Machol,
  Omer, Lander, Aiden, "A 3D map of the human genome at kilobase resolution
  reveals principles of chromatin looping", *Cell* 159(7):1665-1680 (2014),
  doi:10.1016/j.cell.2014.11.021. GM12878 in situ Hi-C, primary + replicate
  combined, reads with MAPQ >= 30 (`_combined_30.hic`, 39.9 GB). The Arrowhead
  list is the paper's own domain call on the same primary+replicate map.
- The suggested `https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic`
  now returns HTTP 403 (AccessDenied), so the GEO copy was used. Only the needed
  byte ranges were read (HTTP Range requests, ~31 MB in total); the whole file
  was not downloaded.
- Arrowhead file: 248,065 bytes, MD5 `42a23da326be4379b042ae9658ed6a0b`,
  9,274 domains genome-wide.
- Terms: GEO is a public NCBI archive. NCBI states it places no restrictions on
  use or distribution of GEO data, but notes submitters may hold patent,
  copyright or other rights in some data (see the NCBI/GEO disclaimer,
  https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). No licence file is
  attached to GSE63525. Cite Rao et al. 2014 and GSE63525 when using these
  graphs.
- Reader: the `.hic` file is parsed by our own pure-Python reader in
  `hic_build.py` (the `hicstraw` package needs libcurl headers and did not build
  on Windows). Integrity check: after KR balancing, every retained row sums to
  the same value (coefficient of variation < 2e-8 on all six chromosomes), which
  only holds if all matrix blocks were read and decoded correctly.

## Construction

1. **Matrix.** Intra-chromosomal raw counts at 25 kb (BP) for the chromosome,
   symmetrised. Normalisation: the file's own KR vector for that chromosome and
   resolution, `obs[i,j] = raw[i,j] / (kr[i] * kr[j])`.
2. **Bins dropped.** Bins whose KR factor is NaN or <= 0, or whose raw row is
   all zero (centromere, p-arm gaps, unmappable regions). The remaining bins
   are renumbered 0..n-1 in genomic order; `nodes.tsv` gives the mapping.
3. **Edges (exact rule).** For each retained bin i, take the k retained bins
   j != i with the largest `obs[i,j]` (only `obs > 0`; ties broken by smaller
   genomic distance, then smaller j). The undirected, unweighted edge set is the
   union of these choices over all i (so degree >= k, not = k). No self-loops,
   no duplicates, each edge written once as `u<TAB>v` with u < v.
4. **Groups.** Each Arrowhead domain `[x1, x2)` on the chromosome becomes the
   set of retained bins whose midpoint `start + 12.5 kb` lies in `[x1, x2)`.
   Groups with fewer than 3 nodes are dropped (7 in total: chr19 4, chr20 1,
   chr21 2); exact duplicate node sets would also be dropped (there were none).
   Nested and partly overlapping domains are kept as given.

Why KR-observed and not observed/expected: the same top-k rule on
observed/expected (`obs / mean obs at that genomic distance`) picks almost only
far-apart pairs: median edge length 400-1700 bins (10-40 Mb), fewer than 0.5%
of edges at <= 5 bins, and only 0.3-3% of edges inside a domain. At long range
the expected count is tiny, so O/E there is dominated by sampling noise and by
compartment (A/B) signal, not domains. Those numbers are in `summary.json`
(`oe_rule_*`) but no O/E graphs were written.

## Reproduce

Needs numpy and networkx (we used a venv: Python 3.12, numpy 2.5.4, networkx 3.7).
Run from a folder outside the repo, with the download in its own folder:

    python -I paper/analysis/hic_build.py fetch --cache <cache_dir> --chroms 17 18 19 20 21 22
    python -I paper/analysis/hic_build.py build --cache <cache_dir>         --domains <dir>/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz         --out paper/analysis/data_hic

`fetch` reads ~31 MB of byte ranges from the GEO `.hic` file (about 1 min).

## Files

Each folder `gm12878_chr<C>_25kb_k<K>/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1.
- `groups.txt` - one ground-truth group (Arrowhead domain) per line, space-separated node ids, sorted.
- `nodes.tsv` - header `id chrom start`; `start` is the bin start in bp (hg19); bin = [start, start + 25000).

`summary.json` holds every number below plus extra diagnostics.

## Counts

Per chromosome (same for all k):

| chrom | bins total | bins kept (n) | Arrowhead domains | groups kept | group size min-median-max | frac nodes in >=1 group | nodes in >1 group | nested pairs | partial-overlap pairs |
|---|---|---|---|---|---|---|---|---|---|
| chr17 | 3252 | 3083 | 411 | 411 | 3-7-80 | 0.757 | 1068 | 282 | 29 |
| chr18 | 3128 | 2991 | 229 | 229 | 3-8-81 | 0.605 | 616 | 130 | 9 |
| chr19 | 2367 | 2222 | 317 | 313 | 3-6-62 | 0.708 | 777 | 220 | 17 |
| chr20 | 2523 | 2382 | 249 | 248 | 3-7-76 | 0.778 | 568 | 138 | 5 |
| chr21 | 1928 | 1384 | 96 | 94 | 3-7-81 | 0.532 | 176 | 41 | 3 |
| chr22 | 2055 | 1353 | 177 | 177 | 3-7-43 | 0.759 | 518 | 130 | 14 |

Per graph. "within" = fraction of edges whose two ends share at least one
group. "band null" = the same top-k rule with score `1/(1+|i-j|)` (pure
genomic distance, no Hi-C data), shown to separate domain signal from the
distance decay. All graphs are connected except three at k = 8: chr18 and
chr19 split into their two arms at the centromere (617 + 2374 and 971 + 1251
nodes), and chr20 into its two arms plus an 11-node pericentromeric island
(starting at 29.40 Mb).

| graph | n | m | mean degree | clustering | frac edges at dist <= 5 bins | groups | within | within, band null |
|---|---|---|---|---|---|---|---|---|
| gm12878_chr17_25kb_k8 | 3083 | 13787 | 8.94 | 0.697 | 0.894 | 411 | 0.640 | 0.602 |
| gm12878_chr17_25kb_k16 | 3083 | 27684 | 17.96 | 0.740 | 0.545 | 411 | 0.519 | 0.493 |
| gm12878_chr17_25kb_k32 | 3083 | 54538 | 35.38 | 0.752 | 0.280 | 411 | 0.335 | 0.339 |
| gm12878_chr18_25kb_k8 | 2991 | 13051 | 8.73 | 0.675 | 0.928 | 229 | 0.536 | 0.500 |
| gm12878_chr18_25kb_k16 | 2991 | 26619 | 17.80 | 0.731 | 0.553 | 229 | 0.452 | 0.424 |
| gm12878_chr18_25kb_k32 | 2991 | 53378 | 35.69 | 0.751 | 0.279 | 229 | 0.312 | 0.310 |
| gm12878_chr19_25kb_k8 | 2222 | 9927 | 8.94 | 0.702 | 0.910 | 313 | 0.583 | 0.549 |
| gm12878_chr19_25kb_k16 | 2222 | 19768 | 17.79 | 0.735 | 0.552 | 313 | 0.467 | 0.445 |
| gm12878_chr19_25kb_k32 | 2222 | 39390 | 35.45 | 0.736 | 0.280 | 313 | 0.301 | 0.306 |
| gm12878_chr20_25kb_k8 | 2382 | 10564 | 8.87 | 0.700 | 0.909 | 248 | 0.665 | 0.629 |
| gm12878_chr20_25kb_k16 | 2382 | 21496 | 18.05 | 0.738 | 0.543 | 248 | 0.554 | 0.528 |
| gm12878_chr20_25kb_k32 | 2382 | 42912 | 36.03 | 0.747 | 0.276 | 248 | 0.392 | 0.389 |
| gm12878_chr21_25kb_k8 | 1384 | 6115 | 8.84 | 0.672 | 0.902 | 94 | 0.460 | 0.425 |
| gm12878_chr21_25kb_k16 | 1384 | 12384 | 17.90 | 0.720 | 0.542 | 94 | 0.374 | 0.353 |
| gm12878_chr21_25kb_k32 | 1384 | 24617 | 35.57 | 0.728 | 0.276 | 94 | 0.248 | 0.253 |
| gm12878_chr22_25kb_k8 | 1353 | 6109 | 9.03 | 0.694 | 0.880 | 177 | 0.638 | 0.606 |
| gm12878_chr22_25kb_k16 | 1353 | 12207 | 18.04 | 0.733 | 0.535 | 177 | 0.542 | 0.508 |
| gm12878_chr22_25kb_k32 | 1353 | 24143 | 35.69 | 0.734 | 0.274 | 177 | 0.374 | 0.370 |

Clustering is networkx `average_clustering` on the unweighted graph.

## Caveats

- The graphs are close to a genomic "band": at k = 8, ~90% of edges join bins
  <= 5 bins (125 kb) apart. Most of the structure is distance decay along the
  chain. Edges fall inside domains only 2-4 points more often than for the
  pure-distance band null at k = 8 and 16, and no more often at k = 32. The
  domain signal is real but weak, and much of any detection score on these
  graphs may come from genomic contiguity alone.
- Domains are small at 25 kb (median 6-8 bins), many sit inside larger ones,
  and 22-47% of nodes are in no domain, so the ground truth is a partial
  cover, not a partition.
- Arrowhead was run on the 5-10 kb map; mapping to 25 kb bins by midpoint
  moves domain edges by up to 12.5 kb.
- One cell line, one normalisation (KR), one resolution, one domain caller.
