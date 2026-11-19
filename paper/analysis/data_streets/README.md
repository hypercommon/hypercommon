# Street graphs of city centres with named streets as ground truth

24 graphs: 12 city centres, each as `<city>_r2` and `<city>_r3`, plus a
distance-only control `<name>_null` for every graph (48 folders in all).
Nodes are points along the streets, edges join points that are close
*along a street*, and the groups are the named streets. Groups are long thin
chains that cross each other and share nodes: a point at a crossing belongs
to every street through it. Built by `paper/analysis/streets_build.py`
(run with `python -I`). No community detection was run.

## Source and licence

- **Data: (c) OpenStreetMap contributors**, available under the Open Database
  Licence (ODbL) 1.0, https://www.openstreetmap.org/copyright and
  https://opendatacommons.org/licenses/odbl/1-0/. These folders are a derived
  database of OpenStreetMap: if you share them, keep this attribution and
  share them under ODbL 1.0.
- Fetched through the public Overpass API, `https://overpass-api.de/api/interpreter`
  (Overpass API 0.7.62.11 87bfad18), no account or login, on 2026-11-15
  (about 21:02-21:19 UTC). One request at a time, 15 s between requests,
  User-Agent `hypercommon-benchmark/1.0 (research graph build; python-urllib)`,
  retry with doubling back-off (30 s, 60 s, ...) on errors. The server
  answered HTTP 502 one or two times for 6 of the 12 cities; every city
  succeeded on a retry. No `remark` (Overpass runtime error) in any answer.
- Each query is pinned to the OSM state `[date:"2026-11-15T00:00:00Z"]`
  (Overpass attic data), so running the same query again should give the
  same ways. `timestamp_osm_base` in the answers is the server's database
  time (2026-11-15T21:02Z-21:18Z), not the data date.
- The raw answers were saved untouched in their own folder (here
  `D:/tmp/streets_dl/`, not in the repository) and only parsed as JSON.

### Queries

One query per city, a square box of side 3700 m around the centre (big
enough for every box the size rule can choose). The Manhattan query, in full:

```
[out:json][timeout:300][date:"2026-11-15T00:00:00Z"];
way["highway"~"^(primary|secondary|tertiary|residential|unclassified|living_street|pedestrian)$"]["name"](40.7382626,-74.0059633,40.7715374,-73.9620367);
out body geom;
```

The other cities differ only in the bounding box (south, west, north, east),
which `streets_build.py` computes from the centre in `CITIES`:

| city | centre (lat, lon) | Overpass box (s, w, n, e) | answer bytes | answer sha256 | query sha256 |
|---|---|---|---|---|---|
| manhattan_midtown | 40.7549, -73.9840 | 40.7382626, -74.0059633, 40.7715374, -73.9620367 | 1,914,531 | `ed4f1f89ab40b8489a6c0d684b514547dbd59f9efe8bf41911c125a4f952aa4a` | `4c16e1924fc50393d2076ceaa29d78999caa05d4955420abbe63e35daae0c31c` |
| chicago_loop | 41.8807, -87.6298 | 41.8640626, -87.6521460, 41.8973374, -87.6074540 | 2,177,606 | `6d13069944ac06d288cf62eb90a7036e6d07fa41b6e491f603bcafd68812282a` | `308341be9a06e6b2f0ab84bce20b7eb12df108cbd3d4b49f3c9541ac3e77e402` |
| barcelona_eixample | 41.3917, 2.1649 | 41.3750626, 2.1427229, 41.4083374, 2.1870771 | 3,721,086 | `6a3bef3fe705f820a45ce8756f362365ec41e07dce2c3361fbfe38b37e81f619` | `2cca22fb02ec0dd7356a3d31761459ad25c89c16fd037d04488b2056c786fe40` |
| paris_centre | 48.8589, 2.3469 | 48.8422626, 2.3216119, 48.8755374, 2.3721881 | 3,532,025 | `363da06360ab5e0bf22264c956aea2199a9e9ded65e366909f98296d4c7ffa21` | `59e221bf8cf7819981a194638f43468f377d6c772eb2e49b38eb3431391f753f` |
| rome_centre | 41.8986, 12.4769 | 41.8819626, 12.4545477, 41.9152374, 12.4992523 | 2,905,304 | `880d0b022f4e44d9d2dd5afdbfded36100b5069979ad36a17bf6e90c3e84fbaa` | `522f585b8560ca070954883a7f356f82f69be7281c4fd3bdb11cd2b0a09039de` |
| tokyo_ginza | 35.6717, 139.7650 | 35.6550626, 139.7445199, 35.6883374, 139.7854801 | 962,276 | `952d4c1e275f3a4a90592f0b861288f58c8f1e2dc18fdad313ed23d07d008427` | `6719abc5ee6c0b3012df609c8223ac77c3ae91ecb8b30219ff52d1e79d6f1342` |
| kyoto_centre | 35.0060, 135.7620 | 34.9893626, 135.7416880, 35.0226374, 135.7823120 | 1,020,195 | `861374be6f96afecdde2a7d7103e56478c06299afd688e9b2cf01de3a361d61a` | `5b019f7ee7fd7a106c74271a8e770cd7d4d8d80300b31e0f42127658038ce8d5` |
| amsterdam_centre | 52.3731, 4.8926 | 52.3564626, 4.8653486, 52.3897374, 4.9198514 | 2,905,604 | `3aa5faeb187b0946930199aa7ab7ed91f45d53cb48086a212a8aef16bc684cf1` | `67a2e064686e8145a9497fdb10e39dc750debd57dc942a8d7269ec390ff69b4e` |
| london_centre | 51.5130, -0.1300 | 51.4963626, -0.1567338, 51.5296374, -0.1032662 | 4,799,590 | `1a085fb6f69f8d6de13bd9c2eb94e42a380c92c56a524061ffdfc135d683dff7` | `2423f4d6ba04778a74fd9674f52894e627c8c98e3ea4b8f9f16c0aed8190416e` |
| prague_old_town | 50.0870, 14.4210 | 50.0703626, 14.3950698, 50.1036374, 14.4469302 | 2,899,001 | `fd854cfdf81e15a77ac1180fd0e7cf8c080d5b99c76de58ec872a9ea3ef3b6ba` | `5c7164ec15c5bed9518bc552e63e8fcb134645cff5543a3da5296f2b5fc059c1` |
| melbourne_cbd | -37.8136, 144.9631 | -37.8302374, 144.9420402, -37.7969626, 144.9841598 | 3,409,703 | `89ddb78042ebdc8b4e180f2f954945bdd7f778af59906058ecbb167800d5a7c3` | `4b3489e2e8b4f09a502b94c9033779775973be726e0fc3fba796088b57cf2fcb` |
| buenos_aires_centro | -34.6037, -58.3816 | -34.6203374, -58.4018131, -34.5870626, -58.3613869 | 1,436,556 | `c71cb617de0125673bbde16ec22216b8a9d8080afbc1fd91f0350aac349949d1` | `33343631b228bcfa0e4e74cc400ac1a5ed8d38884aba5351a623c81815f2db13` |

Query files are saved as `<city>.overpassql` (LF line ends, with the
trailing newline) and answers as `<city>.json`, side by side.

## Construction

All rules below were fixed in `streets_build.py` before any graph was built
or any statistic looked at. Nothing was tuned afterwards.

1. **Areas.** 12 fixed centres (table above): grid cities (Manhattan,
   Chicago, Barcelona Eixample, Kyoto, Melbourne, Buenos Aires) and irregular
   old towns (Paris, Rome, Amsterdam, London, Prague; Tokyo Ginza is a grid
   with many unnamed streets). Boxes are squares in a local equirectangular
   projection around the centre: `x = R cos(lat0) dlon`, `y = R dlat`,
   R = 6,371,008.8 m. This projection is used for all lengths and distances.
2. **Box size rule.** Sides `s_j = 1500 m * 1.25^j`, j in -4..4. Start at
   j = 0; if n (step 5) > 5000 go to j - 1, if n < 2000 go to j + 1; stop at
   the first j with 2000 <= n <= 5000. (Leaving -4..4 would drop the city;
   jumping over the range would take the j closest to it in log n. Neither
   happened.) Result: 1500 m for 10 cities; Kyoto (n = 1628 at 1500 m) and
   Buenos Aires (n = 1871) went to 1875 m.
3. **Ways.** `highway` in {primary, secondary, tertiary, residential,
   unclassified, living_street, pedestrian} with a `name` tag. Ways tagged
   `area=yes` (pedestrian squares drawn as outlines, not streets) are left
   out: 1-189 per city (189 in Barcelona, which maps many squares and
   street sections as areas).
4. **Clip.** A way segment (two consecutive OSM nodes) is kept iff both end
   nodes are inside the box. Segments crossing the border are dropped, not
   cut. Two ways using the same node pair give one segment (never happened).
5. **Nodes.** The OSM nodes that end a kept segment, plus interior points:
   a segment of length L gets `ceil(L / 30 m) - 1` evenly spaced interior
   points, so neighbouring points are at most 30 m apart. A shared OSM node
   (a crossing) stays one node. Positions are linear in lat/lon. Node ids:
   OSM nodes in OSM-id order, then interior points in order of (smaller OSM
   id, larger OSM id, index from the smaller end).
6. **G1.** Consecutive points along each kept segment are linked (a union of
   paths; degree 2 inside a street, 1 at a dead end, 3-4+ at a crossing).
7. **Graphs.** `<city>_r2` = G1 squared (all pairs within 2 hops of G1),
   `<city>_r3` = G1 cubed (within 3 hops). No self-loops. On a plain stretch
   of street this gives degree 4 (r2) or 6 (r3) and triangles, so each
   street becomes a thick chain; at a crossing the hops reach into the other
   street.
8. **Control `<name>_null`.** Same nodes, same `groups.txt` and `nodes.tsv`.
   k = floor(median degree of the real graph) (4 for every r2; 7 or 8 for
   r3). Each node keeps its k nearest other nodes by straight-line distance,
   ties by smaller id; union over nodes, so mean degree is a bit above k.
   No street information, only geometry.
9. **Groups.** All kept segments of ways with the same exact `name` string
   are merged; the group is the set of their points (end nodes and interior
   points). Groups with < 3 points are dropped; exact duplicate point sets
   would be kept once (none). Written in code-point order of the name; the
   names themselves are not written.

Why this rule: the groups are thin, long and cross, unlike the blob-shaped
groups of the other real families. Squaring/cubing G1 keeps every edge on a
street (or across one crossing) but gives each chain triangles; the k-nearest
control has the same nodes and a similar degree but links points across
blocks, ignoring the streets.

## Reproduce

Needs numpy, scipy and networkx (we used a separate venv: Python 3.12.10,
numpy 2.5.4, scipy 1.18.1, networkx 3.7). From `paper/`:

    python -I analysis/streets_build.py fetch --dl <own empty folder>
    python -I analysis/streets_build.py build --dl <same folder> --out analysis/data_streets

`fetch` skips a city whose `.json` and identical `.overpassql` are already
there. `build` takes about 75 s; a second build from the same answers gave
byte-identical folders. To check you have the same input, compare the
sha256 values above (a fresh fetch may differ in bytes only if the server
answers attic queries differently; the ways should be the same).

## Files

Each folder `<city>_r<R>/` and `<city>_r<R>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v.
- `groups.txt` - one street per line, space-separated sorted node ids, size >= 3.
- `nodes.tsv` - header `id lat lon` (WGS84, 7 decimals). No names, no OSM ids.

`summary.json` has one row per real graph, with the control's numbers under
`null_*`. Fields: n, m, mean / min / median / max degree, clustering
(networkx `average_clustering`), n_components, largest_component,
median / max edge length (m), `frac_edges_within_group` (fraction of edges
whose ends share at least one street), side_m, box, n_osm_nodes,
n_interior_points, street_length_km, the same G1 numbers under `g1_*`,
n_groups and drop counts, group size min / median / mean / max,
`frac_nodes_in_group`, `groups_per_node_mean` (memberships per node),
`frac_nodes_in_ge2_groups`, `groups_per_node_max`, null_k, the size-rule
walk, the answer sha256 and `osm_base`. `cities.json` has one row per city
(centre, size-rule walk, named-way counts, download metadata).

## Counts

| graph | box side (m) | n | m | mean deg | clustering | comps | groups | size min-med-max | memb/node | frac >=2 groups | within | null k | null m | within, null |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| manhattan_midtown_r2 | 1500 | 2042 | 4791 | 4.69 | 0.522 | 7 | 54 | 4-37-226 | 1.079 | 0.073 | 0.883 | 4 | 5297 | 0.754 |
| manhattan_midtown_r3 | 1500 | 2042 | 7912 | 7.75 | 0.600 | 7 | 54 | 4-37-226 | 1.079 | 0.073 | 0.795 | 8 | 9891 | 0.614 |
| chicago_loop_r2 | 1500 | 2253 | 5394 | 4.79 | 0.515 | 7 | 68 | 4-34.5-102 | 1.100 | 0.084 | 0.900 | 4 | 5652 | 0.797 |
| chicago_loop_r3 | 1500 | 2253 | 9131 | 8.11 | 0.592 | 7 | 68 | 4-34.5-102 | 1.100 | 0.084 | 0.810 | 8 | 11110 | 0.683 |
| barcelona_eixample_r2 | 1500 | 2827 | 6794 | 4.81 | 0.520 | 13 | 87 | 3-15-206 | 1.098 | 0.093 | 0.875 | 4 | 7239 | 0.742 |
| barcelona_eixample_r3 | 1500 | 2827 | 11404 | 8.07 | 0.594 | 13 | 87 | 3-15-206 | 1.098 | 0.093 | 0.781 | 8 | 14063 | 0.598 |
| paris_centre_r2 | 1500 | 3270 | 7801 | 4.77 | 0.510 | 5 | 248 | 3-10-113 | 1.140 | 0.125 | 0.858 | 4 | 8369 | 0.701 |
| paris_centre_r3 | 1500 | 3270 | 12958 | 7.93 | 0.584 | 5 | 248 | 3-10-113 | 1.140 | 0.125 | 0.748 | 8 | 15918 | 0.558 |
| rome_centre_r2 | 1500 | 3608 | 8190 | 4.54 | 0.521 | 25 | 343 | 3-8-89 | 1.131 | 0.117 | 0.883 | 4 | 9029 | 0.789 |
| rome_centre_r3 | 1500 | 3608 | 13379 | 7.42 | 0.602 | 25 | 343 | 3-8-89 | 1.131 | 0.117 | 0.785 | 7 | 15626 | 0.688 |
| tokyo_ginza_r2 | 1500 | 2085 | 4897 | 4.70 | 0.521 | 15 | 82 | 3-16.5-275 | 1.098 | 0.095 | 0.904 | 4 | 5119 | 0.816 |
| tokyo_ginza_r3 | 1500 | 2085 | 8195 | 7.86 | 0.600 | 15 | 82 | 3-16.5-275 | 1.098 | 0.095 | 0.822 | 7 | 8763 | 0.732 |
| kyoto_centre_r2 | 1875 | 2948 | 7328 | 4.97 | 0.503 | 4 | 64 | 3-40.5-195 | 1.114 | 0.115 | 0.870 | 4 | 7306 | 0.790 |
| kyoto_centre_r3 | 1875 | 2948 | 12513 | 8.49 | 0.571 | 4 | 64 | 3-40.5-195 | 1.114 | 0.115 | 0.764 | 8 | 13936 | 0.631 |
| amsterdam_centre_r2 | 1500 | 3375 | 7927 | 4.70 | 0.517 | 15 | 268 | 3-6-182 | 1.138 | 0.128 | 0.860 | 4 | 8373 | 0.764 |
| amsterdam_centre_r3 | 1500 | 3375 | 13137 | 7.78 | 0.591 | 15 | 268 | 3-6-182 | 1.138 | 0.128 | 0.752 | 8 | 16370 | 0.642 |
| london_centre_r2 | 1500 | 4878 | 11361 | 4.66 | 0.512 | 9 | 333 | 3-12-157 | 1.123 | 0.113 | 0.886 | 4 | 12476 | 0.758 |
| london_centre_r3 | 1500 | 4878 | 18723 | 7.68 | 0.590 | 9 | 333 | 3-12-157 | 1.123 | 0.113 | 0.790 | 8 | 23925 | 0.614 |
| prague_old_town_r2 | 1500 | 2709 | 6151 | 4.54 | 0.511 | 7 | 165 | 3-15-101 | 1.107 | 0.093 | 0.907 | 4 | 6779 | 0.801 |
| prague_old_town_r3 | 1500 | 2709 | 9971 | 7.36 | 0.594 | 7 | 165 | 3-15-101 | 1.107 | 0.093 | 0.826 | 7 | 11689 | 0.695 |
| melbourne_cbd_r2 | 1500 | 3574 | 8378 | 4.69 | 0.517 | 25 | 117 | 3-8-329 | 1.092 | 0.090 | 0.912 | 4 | 8963 | 0.817 |
| melbourne_cbd_r3 | 1500 | 3574 | 14055 | 7.87 | 0.595 | 25 | 117 | 3-8-329 | 1.092 | 0.090 | 0.829 | 7 | 15637 | 0.762 |
| buenos_aires_centro_r2 | 1875 | 2817 | 6926 | 4.92 | 0.509 | 3 | 61 | 5-57-118 | 1.097 | 0.090 | 0.864 | 4 | 7323 | 0.761 |
| buenos_aires_centro_r3 | 1875 | 2817 | 11671 | 8.29 | 0.581 | 3 | 61 | 5-57-118 | 1.097 | 0.090 | 0.762 | 8 | 13516 | 0.611 |

"comps" = connected components of the real graph (same as of G1). "within"
= fraction of edges whose two ends share a street.

Summary over the 12 cities: memberships per node 1.08-1.14, and 7-13% of
nodes are in >= 2 streets. Within-street fraction, real vs control:
r2 0.858-0.912 (median 0.883) vs 0.701-0.817 (median 0.776), gap 8-16
points; r3 0.748-0.829 (median 0.788) vs 0.558-0.762 (median 0.636), gap
7-19 points. The real graph is above its control in every city.

## Caveats

- **Overlap is real but thin.** Only crossing nodes are shared, and
  densifying adds many single-street points, so memberships per node are
  only 1.08-1.14 and a node is in at most 3-5 streets. Overlap sits at
  points, not along stretches.
- **Not connected.** The real graphs have 3-25 components: named streets
  linked only through unnamed ways, footways, or ways cut at the box edge
  fall apart. Amsterdam is the worst: the largest component holds only 1864
  of 3375 nodes (55%), probably because many bridges and canal-side links
  are unnamed or not in the highway list. Elsewhere the largest component
  holds 91-99% of the nodes. Controls have 1-9 components.
- **Groups are not always connected.** 78-94% of the groups per city are
  connected in G1; the rest are a street name split into
  pieces (cut by the box, a square, a section with a different name, or two
  same-named streets). Divided roads mapped as two one-way ways with one
  name give one group made of two parallel chains.
- **Group sizes are skewed**: 3 to 329 points, medians 6-57. Old towns
  (Rome, Amsterdam, Paris) have many short streets; grid cities have fewer,
  longer ones.
- **Degree is low and fixed by the rule**: mean 4.5-5.0 (r2) and 7.4-8.5
  (r3). Some nodes have degree 1 (two-point fragments).
- **OSM ways are not split at the same density everywhere**: OSM nodes make
  up 52-86% of the points (curvy old-town streets have many shape nodes, a
  grid has few), so point spacing varies (median r2 edge length 13.6-27.1 m,
  `edge_length_median_m`).
- **Tokyo**: most small Japanese streets have no `name` tag, so the Ginza
  graph is only the named avenues and lanes (n = 2085 at 1500 m, near the
  lower bound).
- **Name matching is exact**: "West 42nd Street" and "East 42nd Street"
  are two groups; spelling variants are not merged.
- OSM data change over time; the date-pinned query and the sha256 values
  above fix the input.
