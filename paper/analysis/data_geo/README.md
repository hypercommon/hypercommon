# Commuting-flow graphs of areal units with regions as ground truth

10 graphs plus 10 distance-only controls (20 folders). Nodes are small areal
units on a map; edges come from commuting flows between them; ground-truth
groups are regions. Built by `paper/analysis/geo_build.py` (run with `python -I`).

| family | nodes | groups | k |
|---|---|---|---|
| `us_county_state_k{8,16}` | 3,109 US counties (48 contiguous states + DC; Connecticut planning regions) | states | 8, 16 |
| `us_county_cz_k{8,16}` | same graph | USDA ERS 2020 commuting zones (see the circularity caveat) | 8, 16 |
| `ew_ltla_region_k{8,16}` | 331 local authorities (LTLA = LAD 2022) of England and Wales | 9 English regions + Wales | 8, 16 |
| `ew_msoa_north_lad_k{8,16}` | 3,279 MSOAs (2021) of North East, North West, Yorkshire and the Humber, East Midlands, West Midlands | local authority districts (LAD 2022) | 8, 16 |
| `ew_msoa_south_lad_k{8,16}` | 3,985 MSOAs (2021) of East of England, London, South East, South West, Wales | local authority districts (LAD 2022) | 8, 16 |
| `<name>_null` | same nodes and groups as `<name>` | same | same k |

The two US graph sets share their edges (and their controls share edges); only
`groups.txt` differs. England and Wales has 7,264 MSOAs, above the size limit,
so it is split into two halves along region borders.

## Sources

All downloaded on 2026-11-13, without an account or login. Each download went
into its own empty folder; the files were only read as data (the `.xlsx` is
read with openpyxl, no macros run).

| what | URL | file, sha256 |
|---|---|---|
| ACS county-to-county commuting flows, 5-year 2016-2020, Table 1 (sorted by residence) | https://www2.census.gov/programs-surveys/demo/tables/metro-micro/2020/commuting-flows-2020/table1.xlsx (page: https://www.census.gov/data/tables/2020/demo/metro-micro/commuting-flows-2020.html) | `table1.xlsx`, 6,082,366 B, `1dbf52bede8af185ec487e193903e19332c9633cbb3da689afcc382fef7317ee` |
| Census 2022 Gazetteer, counties (internal-point lat/lon) | https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2022_Gazetteer/2022_Gaz_counties_national.zip | zip `7315c7d599e96b220b34402210940032fa4e0278af2b90d8248df1ebbb6405ce`; txt `c33a72d207c78cf7493bbe63d151eb88a8d1235c297b54b59edccdf10af2e8c3` |
| USDA ERS 2020 Commuting Zones (final, updated 9 Jan 2026) | https://www.ers.usda.gov/media/6968/2020-commuting-zones.csv (page: https://www.ers.usda.gov/data-products/commuting-zones-and-labor-market-areas) | `2020-commuting-zones.csv`, 206,760 B, `595e06fea42335504d828cfe1de73044e6df3c8b896ccda777b30e4223f35a6f` |
| ONS Census 2021 origin-destination, ODWP01EW "location of usual residence and place of work" (Nomis bulk) | https://www.nomisweb.co.uk/output/census/2021/odwp01ew.zip (page: https://www.nomisweb.co.uk/sources/census_2021_od) | zip 80,621,150 B `9e32ababfd9f77e353411d399e463f812942440df115a2d6c296c74dbeea70d7`; `ODWP01EW_LTLA.csv` `7bbbf58494cb1a1bc69e00ebacaa281647688bddc84cbce2d28a850c069c1033`; `ODWP01EW_MSOA.csv` `8af475023e18227fdcee3ac4a547d6549b1fcb87138bd93d177a7a698d1a10dd` |
| ONS lookup MSOA (2021) -> BUA (2022) -> LAD (2022) -> region (2022), E&W | https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/MSOA21_BUA22_LAD22_RGN22_EW_LU_v2/FeatureServer/0 (ArcGIS REST `query`, 8 pages of 1,000 rows, attributes only) | 7,264 rows, one per MSOA |
| ONS MSOA (Dec 2021) boundaries, BSC V3, attributes `LAT`, `LONG`, `BNG_E`, `BNG_N` | https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/MSOA_2021_EW_BSC_V3_RUC/FeatureServer/0 (REST `query`, 4 pages, no geometry) | 7,264 rows |

Notes on the sources:
- The 2016-2020 ACS flow table already uses the nine Connecticut planning
  regions (09110-09190) as county equivalents, and so does the ERS 2020 CZ file.
  The 2020 Gazetteer still has the eight old Connecticut counties, so the 2022
  Gazetteer (the first with planning regions) is used for the US points. All
  3,109 contiguous-US + DC units match exactly across the three files.
- The ONS population-weighted MSOA centroid file (hub download) returned HTTP
  403, and its REST service needs a token, so the boundary layer's own `LAT`,
  `LONG` attribute (ONS "centroid or point inside the polygon") is used. That is
  the same kind of point as the US Gazetteer "internal point".
- The Nomis zip download dropped twice mid-file (server does not support
  resume); the third attempt was complete and passes the zip CRC check.

### Licence and terms

- **US Census Bureau** (ACS flows, Gazetteer): works of the US federal
  government, not under copyright in the US (17 U.S.C. 105). Cite "U.S. Census
  Bureau, 2016-2020 American Community Survey 5-year commuting flows" and the
  2022 Gazetteer files.
- **USDA ERS** commuting zones: US federal government work, no licence or use
  terms stated. Requested citation: "U.S. Department of Agriculture, Economic
  Research Service. 2020 Commuting Zones, January 2026." Method: Fowler (2024),
  "New Commuting Zone delineation for the U.S. based on 2020 data",
  *Scientific Data*, with an added contiguity constraint.
- **ONS** census tables and geography: Open Government Licence v3.0
  (https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
  Required attribution: "Source: Office for National Statistics licensed under
  the Open Government Licence v3.0. Contains OS data (c) Crown copyright and
  database right 2026." (the geography files contain both OS and ONS
  intellectual property).

## Construction

### Nodes and points

- **US.** Every 2022 Gazetteer county except Alaska, Hawaii and Puerto Rico,
  sorted by 5-digit FIPS, ids 0..3108. Point = Gazetteer `INTPTLAT`, `INTPTLONG`.
- **E&W LTLA.** The 331 LAD 2022 codes (identical to the census LTLA codes),
  sorted by code. Point = mean of the LTLA's MSOA points weighted by the MSOA's
  resident workers (all ODWP01EW rows of that MSOA, every place-of-work
  category).
- **E&W MSOA.** The MSOAs of the five regions of each half, sorted by code.
  Point = ONS boundary `LAT`, `LONG`.

### Edge rule (exact; fixed before any graph was built)

For one family with node set V (flows to units outside V are ignored):

1. `F[a,b]` = workers who live in unit a and work in unit b, a != b, both in V.
   Self-flow (live and work in the same unit, including working from home) is
   ignored; workplaces abroad, offshore, in Scotland or Northern Ireland, or
   outside the half are dropped.
2. Symmetrise: `S = F + F^T`.
3. Coverage normalisation: `R[a] = sum_b S[a,b]` (the unit's total
   cross-border commuting) and `score[a,b] = S[a,b] / sqrt(R[a] R[b])`.
4. Each unit a keeps the k units b != a with the largest score (only
   `score > 0`; ties by smaller great-circle distance between the points, then
   smaller id). The undirected, unweighted edge set is the union over a, so
   degree is usually >= k but can be lower for a unit with fewer than k
   positive-flow partners (30 US counties at k = 8, 254 at k = 16, none in E&W).
   Each edge is written once as `u<TAB>v`, u < v.

Why this normalisation. It is the commuting analogue of the vanilla-coverage
square-root ("VC_SQRT") balancing used for Hi-C maps: a unit's total
cross-border flow plays the part of a bin's coverage. Ranking partners by raw
flow sends every unit's top-k to the nearest big employment centre, so a few
metro units collect hundreds of edges. Dividing by `sqrt(R[a] R[b])` keeps a
link strong only when it is large relative to both ends' total commuting. Only
the partner's factor changes a unit's own ranking (its own factor is constant
along its row), so this is a mild correction for partner size. Self-flow is
left out of R because it was ignored in step 1. We chose this over the
Tolbert-Sizer measure `S / min(W_a, W_b)` because that measure is the
dissimilarity ERS used to cluster the commuting zones, and over raw flow
because of the hubs. The raw-flow rule was computed only as a diagnostic
(table below); no raw-flow graphs were written.

### Control (`<name>_null`)

Same nodes, `groups.txt` and `nodes.tsv`. Each unit keeps its k nearest units
by great-circle distance between the points (haversine, R = 6371.0088 km; ties
by smaller id), union over units. Pure geometry, no flow signal.

### Groups

One line per group, space-separated node ids, sorted; groups with fewer than 3
nodes are dropped (their nodes are in no group).

- `us_county_state`: state FIPS (first two digits). 49 labels; DC (one unit) is
  dropped, leaving 48 groups. Partition of all other units.
- `us_county_cz`: ERS `CZ2020`, restricted to contiguous-US units. 564 labels;
  82 CZs with 1-2 counties are dropped, leaving 482 groups covering 95.7% of the
  counties. CZs can cross state lines.
- `ew_ltla_region`: `RGN22CD` (Wales is one unit, `W92000004`). 10 groups,
  partition.
- `ew_msoa_*_lad`: `LAD22CD`. North: 137 groups, partition. South: 194 labels;
  City of London and Isles of Scilly (one MSOA each) are dropped, leaving 192.

Group order: by state FIPS, CZ code (numeric), region code or LAD code.

## Files

Each folder `<name>/` and `<name>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, no header.
- `groups.txt` - one ground-truth group per line.
- `nodes.tsv` (UTF-8, header row):
  - US: `id fips name state lat lon cz2020 cz2020_name`
  - LTLA: `id lad22cd name rgn22cd rgn22nm lat lon n_msoa`
  - MSOA: `id msoa21cd name lad22cd lad22nm rgn22nm lat lon bng_e bng_n`

`summary.json` has one row per graph with everything below (keys
`null_*` are the control, `rawflow_rule_*` the raw-flow diagnostic) plus flow
totals.

## Counts

"within" = fraction of edges whose two ends share at least one group.
Clustering is networkx `average_clustering`. Every control is connected.

| graph | n | m | mean degree | max degree | clustering | components | groups | group size min-median-max | frac nodes in group | within | within, control | control clustering | median edge km (control) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| us_county_state_k8 | 3109 | 14484 | 9.32 | 26 | 0.412 | 1 | 48 | 3-64-254 | 1.000 | 0.845 | 0.831 | 0.486 | 56.8 (53.1) |
| us_county_cz_k8 | 3109 | 14484 | 9.32 | 26 | 0.412 | 1 | 482 | 3-6-17 | 0.957 | 0.422 | 0.426 | 0.486 | 56.8 (53.1) |
| us_county_state_k16 | 3109 | 28775 | 18.51 | 85 | 0.355 | 1 | 48 | 3-64-254 | 1.000 | 0.755 | 0.775 | 0.553 | 82.7 (72.9) |
| us_county_cz_k16 | 3109 | 28775 | 18.51 | 85 | 0.355 | 1 | 482 | 3-6-17 | 0.957 | 0.261 | 0.294 | 0.553 | 82.7 (72.9) |
| ew_ltla_region_k8 | 331 | 1525 | 9.21 | 20 | 0.555 | 1 | 10 | 12-31.5-64 | 1.000 | 0.846 | 0.819 | 0.562 | 23.2 (23.0) |
| ew_ltla_region_k16 | 331 | 3117 | 18.83 | 42 | 0.606 | 1 | 10 | 12-31.5-64 | 1.000 | 0.761 | 0.728 | 0.619 | 33.0 (32.2) |
| ew_msoa_north_lad_k8 | 3279 | 16781 | 10.24 | 66 | 0.571 | 1 | 137 | 5-16-132 | 1.000 | 0.853 | 0.773 | 0.550 | 3.6 (2.9) |
| ew_msoa_north_lad_k16 | 3279 | 32216 | 19.65 | 104 | 0.610 | 1 | 137 | 5-16-132 | 1.000 | 0.758 | 0.678 | 0.609 | 4.8 (4.1) |
| ew_msoa_south_lad_k8 | 3985 | 20051 | 10.06 | 42 | 0.574 | 2 | 192 | 7-18-73 | 0.999 | 0.829 | 0.747 | 0.555 | 3.6 (3.1) |
| ew_msoa_south_lad_k16 | 3985 | 38353 | 19.25 | 112 | 0.601 | 2 | 192 | 7-18-73 | 0.999 | 0.715 | 0.629 | 0.612 | 5.0 (4.5) |

Median "within" over the 10 graphs: 0.759 (flow) vs 0.737 (control).

The two components of the south MSOA graphs: the 18 Isle of Wight MSOAs keep
all their top-k partners on the island.

Controls and diagnostics:

| graph | control m | control mean degree | raw-flow rule: max degree | raw-flow rule: within | flow edges also in control |
|---|---|---|---|---|---|
| us_county_state_k8 | 13830 | 8.90 | 92 | 0.823 | 0.741 |
| us_county_cz_k8 | 13830 | 8.90 | 92 | 0.387 | 0.741 |
| us_county_state_k16 | 27539 | 17.72 | 255 | 0.722 | 0.652 |
| us_county_cz_k16 | 27539 | 17.72 | 255 | 0.232 | 0.652 |
| ew_ltla_region_k8 | 1653 | 9.99 | 33 | 0.840 | 0.885 |
| ew_ltla_region_k16 | 3360 | 20.30 | 97 | 0.723 | 0.858 |
| ew_msoa_north_lad_k8 | 16192 | 9.88 | 127 | 0.808 | 0.668 |
| ew_msoa_north_lad_k16 | 31921 | 19.47 | 191 | 0.709 | 0.725 |
| ew_msoa_south_lad_k8 | 19575 | 9.82 | 652 | 0.735 | 0.689 |
| ew_msoa_south_lad_k16 | 38476 | 19.31 | 879 | 0.631 | 0.736 |

## Reproduce

Needs numpy, pandas, openpyxl and networkx (we used a venv: Python 3.12,
numpy 2.5.4, pandas 3.1.1, openpyxl 3.1.5, networkx 3.7). Put the downloads
under one folder `<raw>` with these sub-folders (zips extracted with Python's
`zipfile` into separate `_x` folders):

    <raw>/acs_flows/table1.xlsx
    <raw>/gazetteer2022_x/2022_Gaz_counties_national.txt
    <raw>/ers_cz/2020-commuting-zones.csv
    <raw>/ons_odwp01ew_x/ODWP01EW_LTLA.csv, ODWP01EW_MSOA.csv
    <raw>/ons_lookup/page_{0,1000,...,7000}.json
    <raw>/ons_msoa_pts/page_{0,2000,4000,6000}.json

The two ONS REST folders are the `query` responses for
`where=1=1&outFields=MSOA21CD,MSOA21NM,LAD22CD,LAD22NM,RGN22CD,RGN22NM&orderByFields=ObjectId&resultOffset=<o>&resultRecordCount=1000&returnGeometry=false&f=json`
and
`where=1=1&outFields=MSOA21CD,MSOA21NM,BNG_E,BNG_N,LAT,LONG&orderByFields=FID&resultOffset=<o>&resultRecordCount=2000&returnGeometry=false&f=json`.
Then, from `paper/`:

    python -I analysis/geo_build.py --raw <raw> --out analysis/data_geo

This writes all 20 folders and `summary.json` (about 3 min). `--only us` or
`--only ew` rebuilds one country's folders without rewriting `summary.json`.

## Caveats

- **Commuting zones are partly circular.** ERS built the 2020 CZs by
  constrained hierarchical clustering of the same 2016-2020 ACS county commuting
  flows used for the edges here (dissimilarity `1 - (F_ab + F_ba) / min(labour
  force)`, average linkage, contiguity constraint). Set (b) is therefore partly
  a clustering of the edge data itself. States, regions and local authorities
  are administrative and were not derived from commuting.
- **Geometry dominates.** 65-89% of flow edges are also edges of the
  distance-only control. For the US graphs the flow edges are no more often
  inside a state or CZ than the control's (within differs by -3.2 to +1.4
  points), and US flow graphs have lower clustering than their controls (long
  links to metro centres). The LTLA graphs beat the control by about 3 points,
  the MSOA graphs by about 8 points. Much of any detection score on these graphs
  may come from spatial contiguity alone.
- **Groups are large and coarse** for states (median 64 counties, up to 254 for
  Texas) and regions (12-64 LTLAs); CZs are small (median 6). The LAD groups on
  MSOAs (median 16-18) are closest in scale to the k used.
- **ACS sampling error.** The flows are survey estimates with margins of error;
  many small county pairs have MOE larger than the estimate, and pairs with no
  sampled commuters are absent (hence the low-degree counties).
- **Census 2021 timing.** The England and Wales census day (21 March 2021) fell
  in a COVID-19 lockdown with many workers furloughed or at home; ONS advises
  caution with the 2021 workplace flows. Small cells are also perturbed by
  ONS's cell-key disclosure control.
- **Cut borders.** Each MSOA half uses only flows inside the half, so units on
  the north-south seam lose their cross-seam partners. Flows to Scotland,
  Northern Ireland and abroad are dropped for E&W, and flows to Alaska, Hawaii,
  Puerto Rico and abroad for the US.
- **Points.** US internal points and ONS boundary points are not
  population-weighted; for large rural counties the point can be far from where
  people live. The LTLA point is worker-weighted from MSOA points.
