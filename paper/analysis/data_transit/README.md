# Public-transport graphs with routes as ground truth

Nodes are stations, edges join stations that are close along the network, and
the groups are the routes (lines). A station served by several routes is in
all of them, so the groups are long thin chains that cross and share nodes
(true overlap). Built by `paper/analysis/transit_build.py` (run with
`python -I`). No community detection was run.

## Construction (fixed before any feed was read or any statistic computed)

1. **Feeds.** One official GTFS zip per city, downloaded with curl, no login,
   no API key, each into its own empty folder; the zip is read in place with
   Python's `zipfile` (never extracted, nothing in it executed).
2. **Modes.** rail-like = `route_type` 0, 1, 2 (tram, metro, rail) and the
   extended types 100-199 (railway), 400-499 (urban rail / metro), 900-999
   (tram). bus-like = 3, 11 (trolleybus), 700-799 (bus), 800 (trolleybus).
   Other types (ferry, cable car, funicular, coach, ...) are not used.
   If the rail-like routes give >= 300 stations, the mode is `rail`.
   Otherwise the mode is `railbus`: all rail-like routes plus the bus-like
   routes clipped to a 6 km x 6 km box around the city's centre coordinate
   (fixed list in the script, written before any feed was read).
3. **Trips.** Every trip of a chosen route in `trips.txt` (all service days
   and periods in the feed; union), stops ordered by `stop_sequence`.
4. **Stations.** A stop with a `parent_station` (followed to the top of the
   chain) becomes that parent, at the parent's lat/lon. Stops without a parent
   are merged when they have the same `stop_name` and lie within 150 m
   (great-circle; single linkage), at the mean position. Parent stations are
   not merged with anything else.
5. **G1 (L-space).** Two stations are linked if they are consecutive on at
   least one trip (a step from a station to itself is skipped); union over
   all trips and routes. In `railbus` mode a bus step is kept only when both
   stations are inside the box (no bridging over outside stops); rail-like
   trips are never clipped. Nodes = stations that end at least one G1 edge.
6. **Graphs.** `<city>_<mode>_r2` = square of G1 (link if within 2 hops);
   `<city>_<mode>_r3` = cube of G1 (within 3 hops). This thickens each chain
   so that it has triangles.
7. **Control `<name>_null`.** Same nodes and groups. k = median degree of the
   real graph, rounded half up. Each node keeps its k nearest other nodes by
   great-circle distance (ties by smaller id); union over nodes (symmetrised).
8. **Groups.** Per `route_id`, the set of stations of its trips that are nodes
   (bus routes in `railbus` mode: only stations that end a kept in-box step).
   Keep sets with >= 3 stations; identical sets (different route_ids with the
   same stations) are kept once, under the first route_id in sorted order.
9. **Node ids.** Stations sorted by key (parent `stop_id`, or `~` + smallest
   `stop_id` of a merged cluster), numbered 0..n-1.

The rule above was written (here and in the script docstring) before any
feed was opened. Nothing was changed after the statistics were seen; no city
was swapped or dropped.

## Sources

12 official agency feeds, all downloaded on 2026-11-15 with curl from the
agency's own open-data address, no login, no account, no API key. Each zip
went into its own new empty folder `D:/tmp/transit_dl/<city>/feed.zip`.
The feeds change often (most are rebuilt daily or weekly); the sha256 below
identifies the exact file used. Feed dates are `feed_start_date`-
`feed_end_date` from `feed_info.txt` when the feed has one.

| city | agency | URL | licence | bytes | sha256 | feed dates |
|---|---|---|---|---|---|---|
| nyc | MTA New York City Transit (subway) | https://rrgtfsfeeds.s3.amazonaws.com/gtfs_subway.zip | MTA data feeds, free to use (MTA developer terms; no key for static GTFS) ([terms](https://www.mta.info/developers)) | 5,606,305 | `acb122715bd90aeacf61f4559722ccf0b7de9c63317200a6dd7480246488f1e6` | 20261101-20270515 |
| helsinki | HSL Helsinki Region Transport | https://infopalvelut.storage.hsldev.com/gtfs/hsl.zip | CC BY 4.0 (HSL open data) ([terms](https://www.hsl.fi/en/opendata)) | 88,353,564 | `6f2f7f131bb2a1bb053e23380d0818b55788772d2f610258beeefeda1b492c5d` | 20261113-20270111 |
| berlin | VBB Verkehrsverbund Berlin-Brandenburg | https://www.vbb.de/vbbgtfs | CC BY 4.0 ([terms](https://unternehmen.vbb.de/digitale-services/datensaetze/)) | 84,358,537 | `af8eba8f3e285e762141e82df1fcafd9cd76fff6d06c1ff5cdeaad1da27a0ee1` | not given |
| paris | Ile-de-France Mobilites (IDFM) | https://eu.ftp.opendatasoft.com/stif/GTFS/IDFM-gtfs.zip | ODbL ([terms](https://data.iledefrance-mobilites.fr/explore/dataset/offre-horaires-tc-gtfs-idfm/)) | 114,035,165 | `f02e0025e75eb32e29e258c623b948a35a4e5979c279d7bfe56e310e86b1196f` | not given |
| prague | PID / ROPID Prague Integrated Transport | https://data.pid.cz/PID_GTFS.zip | CC BY 4.0 ([terms](https://pid.cz/o-systemu/opendata/)) | 46,347,739 | `dfeb14923baf0c17c42273267ee8d4a5bc9ae306ba63fde728b71a1d9bde2c5f` | 20261115-20261128 |
| budapest | BKK Centre for Budapest Transport | https://go.bkk.hu/api/static/v1/public-gtfs/budapest_gtfs.zip | public timetable data, no formal licence found (BKK developer page) ([terms](https://bkk.hu/bkk-partnerek/egyeb/gtfs-menetrendi-adatok-programozoknak/)) | 63,684,538 | `b611118e5f9f21bd3463e8ac23ee936cf7058edbcb566b498177969761c65964` | 20261114-20270101 |
| boston | MBTA | https://cdn.mbta.com/MBTA_GTFS.zip | MassDOT Developers License Agreement ([terms](https://www.mbta.com/developers/gtfs)) | 22,280,410 | `3c90dac02874275cd232c9502daa7c62749e81b6e03c4d7756f354df873a4c85` | 20261106-20261212 |
| chicago | Chicago Transit Authority (CTA) | https://www.transitchicago.com/downloads/sch_data/google_transit.zip | CTA Developer License Agreement and Terms of Use ([terms](https://www.transitchicago.com/developers/terms/)) | 51,360,701 | `ddbf306e6a57e419437fe01c4e865372f031c0a9626b93214b418457c21b72b4` | not given |
| rome | Roma Servizi per la Mobilita | https://romamobilita.it/sites/default/files/rome_static_gtfs.zip | CC BY 3.0 IT per the agency's own statements; not restated on the current download page ([terms](https://romamobilita.it/it/tecnologie/open-data)) | 46,305,383 | `5ad900cfd383ed7000ed43a1cb79e42a335d958da2fbcb44f27c77324385f9b4` | not given |
| portland | TriMet | https://developer.trimet.org/schedule/gtfs.zip | TriMet Terms of Use (limited licence to use, reproduce, redistribute the Data) ([terms](https://developer.trimet.org/terms_of_use.shtml)) | 43,467,097 | `e642a21ab45100a77824c5f6c2960dc37b632b8191feefb9dcfc2d94fe4d662d` | 20261011-20270227 |
| oslo | Ruter (via Entur national data) | https://storage.googleapis.com/marduk-production/outbound/gtfs/rb_rut-aggregated-gtfs.zip | NLOD 2.0 ([terms](https://developer.entur.org/stops-and-timetable-data)) | 49,690,748 | `45e04b598daae4d3088630b8365af4632a2d4c7e7913c0bd1d556721f58a179d` | not given |
| seattle | King County Metro | https://metro.kingcounty.gov/GTFS/google_transit.zip | King County Transit Data Terms of Use (attribution: 'Transit scheduling, geographic, and real-time data provided by permission of King County') ([terms](https://kingcounty.gov/en/dept/metro/rider-tools/mobile-and-web-apps)) | 10,795,934 | `ff9b26cf8dc49558b526ad804dc1e41b67007fe84914c60bb3a9779c2e5fa163` | 20261029-20270326 |

Licence notes (checked on the agencies' pages on 2026-11-15):

- **CC BY 4.0**: Helsinki (HSL; cite "HSL" and the delivery date), Berlin
  (VBB), Prague (PID). **NLOD 2.0**: Oslo (Entur). **ODbL**: Paris (IDFM;
  derived databases must stay under ODbL and credit IDFM).
- **Agency terms, free use**: NYC (MTA: "Our data feeds are free to use";
  MTA logos/maps need a licence), Boston (MassDOT Developers License
  Agreement), Chicago (CTA Developer License Agreement, a copy is inside the
  zip), Portland (TriMet terms: limited, revocable licence to use, reproduce
  and redistribute the Data), Seattle (King County terms: limited, revocable
  licence to use, reproduce and redistribute; required attribution "Transit
  scheduling, geographic, and real-time data provided by permission of King
  County").
- **Unclear**: Budapest (BKK publishes the feed openly and calls it public
  data, but its developer page could not be fetched and no formal licence
  was found); Rome (the current open-data page names no licence and reserves
  rights against misuse; older agency material says CC BY 3.0 IT). Check
  these two before redistributing the derived graphs.
- Not used: Montreal STM (download returned 403), Vienna (URL now serves an
  HTML page), Melbourne (310 MB nested zip), feeds that need a key (Sydney,
  Stockholm, Washington, Madrid EMT, ...).

## Counts

All rules as above; no detection was run. "within" = fraction of edges
whose two ends share at least one group; "frac >=2" = fraction of nodes in
two or more groups. The control has the same nodes and groups.

| graph | n | m | mean deg | clustering | comp. | groups | size min-med-max | groups/node | frac >=2 | within | null k | null m | null within | null comp. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| nyc_rail_r2 | 496 | 1366 | 5.51 | 0.588 | 6 | 28 | 4-36.5-71 | 2.06 | 0.647 | 0.964 | 5 | 1481 | 0.467 | 2 |
| nyc_rail_r3 | 496 | 2401 | 9.68 | 0.673 | 6 | 28 | 4-36.5-71 | 2.06 | 0.647 | 0.921 | 9 | 2647 | 0.407 | 2 |
| helsinki_railbus_r2 | 346 | 1196 | 6.91 | 0.571 | 4 | 86 | 3-13-55 | 3.91 | 0.734 | 0.841 | 6 | 1285 | 0.534 | 1 |
| helsinki_railbus_r3 | 346 | 2327 | 13.45 | 0.641 | 4 | 86 | 3-13-55 | 3.91 | 0.734 | 0.708 | 11 | 2378 | 0.452 | 1 |
| berlin_rail_r2 | 1549 | 5483 | 7.08 | 0.579 | 7 | 172 | 3-20.5-82 | 2.73 | 0.530 | 0.810 | 5 | 4777 | 0.718 | 1 |
| berlin_rail_r3 | 1549 | 12612 | 16.28 | 0.642 | 7 | 172 | 3-20.5-82 | 2.73 | 0.530 | 0.595 | 10 | 9662 | 0.602 | 1 |
| paris_rail_r2 | 990 | 4156 | 8.40 | 0.579 | 1 | 57 | 3-21-75 | 1.28 | 0.194 | 0.706 | 6 | 3607 | 0.601 | 1 |
| paris_rail_r3 | 990 | 10541 | 21.29 | 0.632 | 1 | 57 | 3-21-75 | 1.28 | 0.194 | 0.460 | 14 | 8468 | 0.471 | 1 |
| prague_rail_r2 | 1141 | 2819 | 4.94 | 0.535 | 3 | 127 | 3-27-151 | 4.37 | 0.599 | 0.947 | 4 | 2767 | 0.737 | 1 |
| prague_rail_r3 | 1141 | 4765 | 8.35 | 0.617 | 3 | 127 | 3-27-151 | 4.37 | 0.599 | 0.890 | 6 | 4112 | 0.688 | 1 |
| budapest_rail_r2 | 447 | 1061 | 4.75 | 0.569 | 6 | 49 | 6-21-47 | 2.30 | 0.593 | 0.903 | 4 | 1070 | 0.650 | 1 |
| budapest_rail_r3 | 447 | 1740 | 7.79 | 0.645 | 6 | 49 | 6-21-47 | 2.30 | 0.593 | 0.822 | 7 | 1864 | 0.601 | 1 |
| boston_railbus_r2 | 751 | 3392 | 9.03 | 0.589 | 1 | 120 | 3-9.5-57 | 2.13 | 0.362 | 0.703 | 6 | 2799 | 0.532 | 1 |
| boston_railbus_r3 | 751 | 8230 | 21.92 | 0.642 | 1 | 120 | 3-9.5-57 | 2.13 | 0.362 | 0.482 | 14 | 6587 | 0.413 | 1 |
| chicago_railbus_r2 | 594 | 1966 | 6.62 | 0.546 | 2 | 48 | 3-27-53 | 2.23 | 0.401 | 0.866 | 5 | 1829 | 0.496 | 1 |
| chicago_railbus_r3 | 594 | 3913 | 13.18 | 0.616 | 2 | 48 | 3-27-53 | 2.23 | 0.401 | 0.732 | 10 | 3597 | 0.425 | 1 |
| rome_railbus_r2 | 795 | 4664 | 11.73 | 0.589 | 3 | 136 | 3-20-57 | 3.82 | 0.688 | 0.770 | 9 | 4259 | 0.494 | 1 |
| rome_railbus_r3 | 795 | 11249 | 28.30 | 0.637 | 3 | 136 | 3-20-57 | 3.82 | 0.688 | 0.577 | 24 | 11177 | 0.381 | 1 |
| portland_railbus_r2 | 599 | 1939 | 6.47 | 0.546 | 3 | 41 | 3-24-54 | 1.75 | 0.432 | 0.851 | 5 | 1849 | 0.532 | 2 |
| portland_railbus_r3 | 599 | 3847 | 12.84 | 0.612 | 3 | 41 | 3-24-54 | 1.75 | 0.432 | 0.703 | 10 | 3638 | 0.462 | 1 |
| oslo_railbus_r2 | 359 | 1095 | 6.10 | 0.574 | 3 | 57 | 3-11-51 | 2.35 | 0.627 | 0.816 | 5 | 1121 | 0.450 | 1 |
| oslo_railbus_r3 | 359 | 2088 | 11.63 | 0.640 | 3 | 57 | 3-11-51 | 2.35 | 0.627 | 0.658 | 9 | 1971 | 0.395 | 1 |
| seattle_railbus_r2 | 428 | 1639 | 7.66 | 0.590 | 1 | 50 | 7-19.5-62 | 2.53 | 0.481 | 0.843 | 6 | 1538 | 0.471 | 1 |
| seattle_railbus_r3 | 428 | 3492 | 16.32 | 0.636 | 1 | 50 | 7-19.5-62 | 2.53 | 0.481 | 0.666 | 13 | 3236 | 0.371 | 1 |

Per city (G1 = the L-space base graph, not written to disk; "merged away" =
parentless stops absorbed by the same-name 150 m rule; "dup" = route_ids
whose station set equals an earlier route's):

| city | mode | rail-like stations | n | nodes in box | stops used | merged away | G1 m | G1 comp. | groups rail / bus | dropped <3 / dup |
|---|---|---|---|---|---|---|---|---|---|---|
| nyc | rail | 496 | 496 | - | 989 | 0 | 578 | 6 | 28 / 0 | 1 / 0 |
| helsinki | railbus | 268 | 346 | 197 | 8277 | 3550 | 445 | 4 | 47 / 39 | 0 / 85 |
| berlin | rail | 1549 | 1549 | - | 3042 | 375 | 1943 | 7 | 172 / 0 | 2 / 6 |
| paris | rail | 990 | 990 | - | 1938 | 0 | 1294 | 1 | 57 / 0 | 0 / 0 |
| prague | rail | 1141 | 1141 | - | 1546 | 342 | 1253 | 3 | 127 / 0 | 0 / 2 |
| budapest | rail | 447 | 447 | - | 901 | 336 | 476 | 6 | 49 / 0 | 0 / 0 |
| boston | railbus | 262 | 751 | 524 | 7601 | 681 | 1001 | 1 | 22 / 98 | 5 / 19 |
| chicago | railbus | 143 | 594 | 479 | 10937 | 4586 | 739 | 2 | 8 / 40 | 0 / 2 |
| rome | railbus | 183 | 795 | 703 | 8333 | 2891 | 1305 | 3 | 9 / 127 | 7 / 12 |
| portland | railbus | 168 | 599 | 535 | 5961 | 1717 | 725 | 3 | 9 / 32 | 0 / 3 |
| oslo | railbus | 183 | 359 | 241 | 7236 | 0 | 423 | 3 | 11 / 46 | 5 / 22 |
| seattle | railbus | 19 | 428 | 428 | 6267 | 1941 | 560 | 1 | 2 / 48 | 1 / 5 |

Every node is in at least one group (frac_nodes_in_group = 1.0 in all but
NYC, 0.996: the 2-station 42 St shuttle is dropped as < 3 stations).

## Reproduce

Needs pandas, numpy, scipy and networkx (we used a separate venv: Python
3.12.10, pandas 3.1.1, numpy 2.5.4, scipy 1.18.1, networkx 3.7). Put each
feed in `<dl>/<city>/feed.zip` (city names as in the table), then from
`paper/`:

    python -I analysis/transit_build.py build --dl <dl> --out analysis/data_transit
    python -I analysis/transit_build.py null  --out analysis/data_transit   # optional: rewrite only *_null/

`build` takes about 8 minutes (single process; Paris and Berlin dominate)
and writes all 48 folders and `summary.json`; `--cities a,b` rebuilds a
subset and merges it into `summary.json`. `null` rebuilds the controls from
the real folders and gives byte-identical files (checked).

## Files

Each folder `<city>_<mode>_r<p>/` and `<city>_<mode>_r<p>_null/` has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v.
- `groups.txt` - one route (station set) per line, space-separated sorted
  node ids, size >= 3.
- `nodes.tsv` - header `id lat lon` (WGS84 degrees, 6 decimals). No route
  or name information.

`summary.json` has one row per real graph, with the control's numbers under
`null_*`: n, m, mean / min / median / max degree, clustering (networkx
`average_clustering`), n_components, `frac_edges_within_group`, `null_k`,
the group counts (`n_groups`, rail / bus, dropped small / duplicate), group
size min / median / mean / max, `groups_per_node_mean`,
`frac_nodes_ge2_groups`, `groups_per_node_max`, `frac_nodes_in_group`, the
G1 numbers (`g1_*`), the station-merging counts, the mode decision
(`n_rail_stations`, `mode`, `n_nodes_in_box`), the route types in the feed,
and the feed's URL, licence, size, sha256 and `feed_info`.

## Caveats

- **Real vs control.** The real graph is above the control on "within" in
  22 of 24 graphs (by 0.07-0.50). The two exceptions are the big regional
  rail networks at r3: Berlin 0.595 vs 0.602 and Paris 0.460 vs 0.471. There
  the cube links across many lines at every interchange, while in the
  suburbs a station's nearest stations are mostly on its own line, so the
  geometric control is already high (0.60-0.74 at r2 for Berlin, Paris,
  Prague). Rail-only networks with long single-line branches make easy
  controls; the r2 graphs keep a clear gap everywhere.
- **Overlap differs a lot.** Mean groups per node is 1.28 in Paris (most
  RER / Transilien stations are on one line) and 3.8-4.4 in Helsinki, Rome
  and Prague. The busiest node is in 41 groups (Rome railbus).
- **Route = route_id.** Some agencies split one line into several route_ids
  (variants, branches, directions) or give one route_id to a line with
  branches. Identical station sets were merged (85 in Helsinki, 22 in Oslo,
  19 in Boston), but near-identical variants (a short-turn or one extra
  stop) stay as separate, almost duplicate groups.
- **Station merging is imperfect.** Interchanges coded as separate parent
  stations are not merged (NYC complexes such as Times Sq - 42 St are
  several parents), so lines that meet there are not linked (transfers are
  not used). This is why NYC G1 has 6 components: IND/BMT (273 stations),
  IRT (154), L (24), 7 (22), Staten Island Railway (21) and the 42 St
  shuttle (2 stations, its group dropped). Parentless stops closer than
  150 m but with different names are not merged either: the bus-heavy
  graphs (Boston, Chicago, Portland, Rome, Seattle) keep 347-643 node pairs
  closer than 150 m (Paris has 1, NYC 50, Berlin 77).
- **Disconnected groups.** In `railbus` mode bus routes are cut at the 6 km
  box; a route that leaves and re-enters it gives a group in pieces, and
  bus routes with < 3 in-box stations are dropped. Overall 0-4.4 % of the
  groups are not connected in the r2 graph (0 % in NYC, Paris, Prague,
  Budapest, Chicago; 1.7 % in Berlin, 4.4 % in Rome).
- **Feed scope.** Berlin (VBB), Prague (PID) and Paris (IDFM) rail-like
  routes cover the whole region, so those graphs are regional, not
  city-only: Berlin spans lat 50.8-54.3, lon 11.0-17.0 (regional trains
  reach into Saxony, Mecklenburg and Poland), Prague lat 49.4-50.6, Paris
  lat 48.0-49.4. Seattle's feed is King County
  Metro only: Link light rail is in Sound Transit's feed, so the rail part is
  just the 2 streetcar lines (19 stations) and the graph is almost all bus.
  Helsinki, Boston, Chicago, Rome, Portland, Oslo and Seattle fall under the
  300-station rail threshold and use `railbus`.
- **Several G1 components.** Pieces not linked to the main network (e.g.
  Staten Island Railway, interchanges split over several parents) are kept;
  the real graphs have 1-7 components, the controls 1-2.
- **Union over the feed period.** All trips of all service days are pooled,
  so temporary diversions and replacement patterns in the feed window add
  edges and stations.
