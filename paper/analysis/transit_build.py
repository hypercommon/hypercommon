"""Build public-transport graphs from official open GTFS feeds: stations are
nodes, routes (lines) are the ground-truth groups. A station served by several
routes is in all of them, so the groups are long thin chains that cross and
share nodes (true overlap).

Source data: one GTFS zip per city, downloaded by hand with curl into its own
empty folder <dl>/<city>/feed.zip (no login, no key, nothing executed). The
zip is read in place with zipfile; it is never extracted. URLs, licences and
sha256 values are in data_transit/README.md and summary.json.

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/transit_build.py build --dl D:/tmp/transit_dl --out analysis/data_transit
  python -I analysis/transit_build.py null  --out analysis/data_transit   # only rewrite *_null/
  (optional: --cities nyc,paris to build a subset; summary.json is then merged)

Rules (fixed before any graph was built or any statistic looked at; no
community detection was run):
  Modes. rail-like = route_type 0, 1, 2 (tram, metro, rail) and the extended
    types 100-199 (railway), 400-499 (urban railway / metro), 900-999 (tram).
    bus-like = route_type 3, 11 (trolleybus), 700-799 (bus), 800 (trolleybus).
    Everything else (ferry, cable car, funicular, coach 200-299, ...) is unused.
    If the rail-like routes give >= MIN_RAIL_N stations (after merging), the
    mode is "rail". Otherwise the mode is "railbus": all rail-like routes plus
    the bus-like routes restricted to a BOX_M x BOX_M box (6 km x 6 km) around
    the city's centre coordinate in CITIES (fixed list, written before any
    download was read).
  Trips. Every trip of a chosen route in trips.txt is used (all service days
    and all feed periods; union). stop_times are ordered by stop_sequence.
  Stations. A stop with parent_station (followed to the top of the chain)
    becomes that parent; the parent's own lat/lon is the station position.
    Stops with no parent are merged when they have the same stop_name and lie
    within MERGE_M = 150 m (great-circle), single linkage; position = mean of
    the merged stops. Parent stations are never merged with each other or with
    parentless stops. The merging is done over the stops used by the chosen
    routes of the mode (for railbus: rail-like plus all bus-like trips, before
    the box clipping).
  G1 (L-space). Stations a, b are linked if they are consecutive on at least
    one trip (consecutive stop_times rows after the station mapping; a step
    from a station to itself is skipped). Union over all trips and routes.
    For a bus-like trip in railbus mode, a step is kept only if both stations
    lie inside the box (no bridging over stops outside it). Rail-like trips
    are never clipped. Nodes = stations that end at least one G1 edge.
  Graphs. <city>_<mode>_r2 = square of G1 (link if distance <= 2 in G1);
    <city>_<mode>_r3 = cube of G1 (distance <= 3).
  Control <name>_null. Same nodes, same groups. k = median degree of the real
    graph rounded half up; each node keeps its k nearest other nodes by
    great-circle distance (ties by smaller id); union (symmetrised).
  Groups. Per route_id, the set of stations of any of its trips that are
    nodes (for bus routes in railbus mode: only the in-box stations, i.e. only
    stations that end a kept step of that route). Keep sets with >= 3
    stations; identical sets (several route_ids with the same stations) are
    kept once (first route_id in sorted order).
  Node ids. Stations sorted by key (parent stop_id, or "~" + the smallest
    stop_id of a merged cluster), then numbered 0..n-1.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import argparse
import hashlib
import io
import json
import math
import shutil
import zipfile
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

MIN_RAIL_N = 300
BOX_M = 6000.0
MERGE_M = 150.0
MIN_GROUP = 3
R_EARTH = 6371008.8
CHUNK = 2_000_000

# Centre coordinate (lat, lon) of each city, fixed before any feed was read:
# the main square / central station named in the comment. Only used when the
# rule falls back to railbus.
CITIES = [
    dict(city="nyc", centre=(40.7128, -74.0060),       # City Hall, Manhattan
         agency="MTA New York City Transit (subway)",
         url="https://rrgtfsfeeds.s3.amazonaws.com/gtfs_subway.zip",
         licence="MTA data feeds, free to use (MTA developer terms; no key for static GTFS)",
         licence_url="https://www.mta.info/developers"),
    dict(city="helsinki", centre=(60.1710, 24.9441),   # Rautatientori
         agency="HSL Helsinki Region Transport",
         url="https://infopalvelut.storage.hsldev.com/gtfs/hsl.zip",
         licence="CC BY 4.0 (HSL open data)",
         licence_url="https://www.hsl.fi/en/opendata"),
    dict(city="berlin", centre=(52.5219, 13.4132),     # Alexanderplatz
         agency="VBB Verkehrsverbund Berlin-Brandenburg",
         url="https://www.vbb.de/vbbgtfs",
         licence="CC BY 4.0",
         licence_url="https://unternehmen.vbb.de/digitale-services/datensaetze/"),
    dict(city="paris", centre=(48.8584, 2.3470),       # Chatelet
         agency="Ile-de-France Mobilites (IDFM)",
         url="https://eu.ftp.opendatasoft.com/stif/GTFS/IDFM-gtfs.zip",
         licence="ODbL",
         licence_url="https://data.iledefrance-mobilites.fr/explore/dataset/offre-horaires-tc-gtfs-idfm/"),
    dict(city="prague", centre=(50.0875, 14.4213),     # Old Town Square
         agency="PID / ROPID Prague Integrated Transport",
         url="https://data.pid.cz/PID_GTFS.zip",
         licence="CC BY 4.0",
         licence_url="https://pid.cz/o-systemu/opendata/"),
    dict(city="budapest", centre=(47.4979, 19.0544),   # Deak Ferenc ter
         agency="BKK Centre for Budapest Transport",
         url="https://go.bkk.hu/api/static/v1/public-gtfs/budapest_gtfs.zip",
         licence="public timetable data, no formal licence found (BKK developer page)",
         licence_url="https://bkk.hu/bkk-partnerek/egyeb/gtfs-menetrendi-adatok-programozoknak/"),
    dict(city="boston", centre=(42.3555, -71.0605),    # Downtown Crossing
         agency="MBTA",
         url="https://cdn.mbta.com/MBTA_GTFS.zip",
         licence="MassDOT Developers License Agreement",
         licence_url="https://www.mbta.com/developers/gtfs"),
    dict(city="chicago", centre=(41.8819, -87.6278),   # State & Madison
         agency="Chicago Transit Authority (CTA)",
         url="https://www.transitchicago.com/downloads/sch_data/google_transit.zip",
         licence="CTA Developer License Agreement and Terms of Use",
         licence_url="https://www.transitchicago.com/developers/terms/"),
    dict(city="rome", centre=(41.8960, 12.4823),       # Piazza Venezia
         agency="Roma Servizi per la Mobilita",
         url="https://romamobilita.it/sites/default/files/rome_static_gtfs.zip",
         licence="CC BY 3.0 IT per the agency's own statements; not restated on the current download page",
         licence_url="https://romamobilita.it/it/tecnologie/open-data"),
    dict(city="portland", centre=(45.5189, -122.6793), # Pioneer Courthouse Sq.
         agency="TriMet",
         url="https://developer.trimet.org/schedule/gtfs.zip",
         licence="TriMet Terms of Use (limited licence to use, reproduce, redistribute the Data)",
         licence_url="https://developer.trimet.org/terms_of_use.shtml"),
    dict(city="oslo", centre=(59.9111, 10.7503),       # Jernbanetorget
         agency="Ruter (via Entur national data)",
         url="https://storage.googleapis.com/marduk-production/outbound/gtfs/rb_rut-aggregated-gtfs.zip",
         licence="NLOD 2.0",
         licence_url="https://developer.entur.org/stops-and-timetable-data"),
    dict(city="seattle", centre=(47.6114, -122.3375),  # Westlake
         agency="King County Metro",
         url="https://metro.kingcounty.gov/GTFS/google_transit.zip",
         licence="King County Transit Data Terms of Use (attribution: 'Transit scheduling, "
                 "geographic, and real-time data provided by permission of King County')",
         licence_url="https://kingcounty.gov/en/dept/metro/rider-tools/mobile-and-web-apps"),
]


def is_rail(t):
    return t in (0, 1, 2) or 100 <= t <= 199 or 400 <= t <= 499 or 900 <= t <= 999


def is_bus(t):
    return t in (3, 11, 800) or 700 <= t <= 799


# ---------------------------------------------------------------- reading
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def member(zf, name):
    """Find a member by basename (some feeds put files in a sub-folder)."""
    hits = [x for x in zf.namelist() if x.rsplit("/", 1)[-1].lower() == name]
    assert len(hits) == 1, (name, hits)
    return hits[0]


def read(zf, name, usecols=None, **kw):
    with zf.open(member(zf, name)) as raw:
        df = pd.read_csv(io.TextIOWrapper(raw, encoding="utf-8-sig"), dtype=str,
                         keep_default_na=False, skipinitialspace=True, **kw)
    df.columns = [c.strip() for c in df.columns]
    if usecols is not None:
        for c in usecols:
            if c not in df.columns:
                df[c] = ""
        df = df[list(usecols)]
    return df


def read_stop_times(zf, trip_code):
    """Rows of stop_times whose trip_id is in trip_code; trip ids -> ints."""
    out = []
    with zf.open(member(zf, "stop_times.txt")) as raw:
        it = pd.read_csv(io.TextIOWrapper(raw, encoding="utf-8-sig"), dtype=str,
                         keep_default_na=False, skipinitialspace=True,
                         chunksize=CHUNK)
        for ch in it:
            ch.columns = [c.strip() for c in ch.columns]
            ch = ch[["trip_id", "stop_id", "stop_sequence"]]
            t = ch["trip_id"].map(trip_code)
            keep = t.notna()
            if keep.any():
                out.append(pd.DataFrame({
                    "t": t[keep].astype(np.int64).to_numpy(),
                    "s": ch.loc[keep, "stop_id"].to_numpy(),
                    "q": pd.to_numeric(ch.loc[keep, "stop_sequence"]).to_numpy()}))
    st = pd.concat(out, ignore_index=True)
    return st.sort_values(["t", "q"], kind="stable").reset_index(drop=True)


def feed_info(zf):
    names = [x.rsplit("/", 1)[-1].lower() for x in zf.namelist()]
    if "feed_info.txt" not in names:
        return {}
    fi = read(zf, "feed_info.txt")
    if not len(fi):
        return {}
    r = fi.iloc[0].to_dict()
    return {k: r[k] for k in ("feed_publisher_name", "feed_version",
                               "feed_start_date", "feed_end_date") if k in r}


# ---------------------------------------------------------------- geometry
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def chord(m):
    return 2 * math.sin(m / (2 * R_EARTH))


def in_box(lat, lon, centre):
    h = BOX_M / 2 / R_EARTH
    dlat = math.degrees(h)
    dlon = math.degrees(h) / math.cos(math.radians(centre[0]))
    return (np.abs(lat - centre[0]) <= dlat) & (np.abs(lon - centre[1]) <= dlon)


# ---------------------------------------------------------------- stations
def station_map(stops, used):
    """stop_id -> station key, and station key -> (lat, lon)."""
    from scipy.spatial import cKDTree
    assert stops["stop_id"].is_unique
    S = stops.set_index("stop_id")
    lat = pd.to_numeric(S["stop_lat"], errors="coerce")
    lon = pd.to_numeric(S["stop_lon"], errors="coerce")
    par = S["parent_station"].str.strip()

    def top(s):
        seen = set()
        while par.get(s, "") and par[s] in S.index and s not in seen:
            seen.add(s)
            s = par[s]
        return s

    smap, pos = {}, {}
    loose = []
    for s in sorted(used):
        assert s in S.index, ("stop_times uses unknown stop", s)
        t = top(s)
        if t != s:
            smap[s] = t
            pos[t] = (float(lat[t]), float(lon[t]))
        else:
            loose.append(s)
    # parentless stops: same name and <= MERGE_M, single linkage
    loose = np.array(loose, dtype=object)
    la = lat[loose].to_numpy(float)
    lo = lon[loose].to_numpy(float)
    nm = S.loc[loose, "stop_name"].str.strip().to_numpy()
    assert np.isfinite(la).all() and np.isfinite(lo).all()
    parent = list(range(len(loose)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    if len(loose):
        tree = cKDTree(xyz(la, lo))
        for i, j in tree.query_pairs(chord(MERGE_M) * 1.001):
            if nm[i] == nm[j] and hav(la[i], lo[i], la[j], lo[j]) <= MERGE_M:
                a, b = find(i), find(j)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    clusters = defaultdict(list)
    for i in range(len(loose)):
        clusters[find(i)].append(i)
    n_merged = 0
    for mem in clusters.values():
        key = "~" + min(loose[i] for i in mem)
        n_merged += len(mem) - 1
        pos[key] = (float(la[mem].mean()), float(lo[mem].mean()))
        for i in mem:
            smap[loose[i]] = key
    return smap, pos, dict(n_stops_used=len(used),
                           n_stops_with_parent=len(used) - len(loose),
                           n_parentless_merged_away=n_merged)


# ---------------------------------------------------------------- G1
def patterns(st, trip_route, smap):
    """Per trip: route and station sequence (consecutive repeats dropped).
    Returns {(route, tuple(seq))} -- trips with the same pattern kept once."""
    stn = st["s"].map(smap).to_numpy()
    t = st["t"].to_numpy()
    cut = np.flatnonzero(np.diff(t)) + 1
    out = set()
    for a, b in zip(np.r_[0, cut], np.r_[cut, len(t)]):
        seq = [stn[a]]
        for x in stn[a + 1:b]:
            if x != seq[-1]:
                seq.append(x)
        out.add((trip_route[t[a]], tuple(seq)))
    return out


def l_space(pats, clip):
    """pats: set of (route, seq, is_bus). clip(station) -> bool for bus steps."""
    E = set()
    route_st = defaultdict(set)
    for r, seq, bus in pats:
        for a, b in zip(seq[:-1], seq[1:]):
            if bus and not (clip(a) and clip(b)):
                continue
            E.add((min(a, b), max(a, b)))
            route_st[r].update((a, b))
    return E, route_st


def powers(n, E1, ks=(2, 3)):
    import scipy.sparse as sp
    e = np.array(E1)
    A = sp.coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    A = ((A + A.T) > 0).astype(np.int8).tocsr()
    out = {}
    R = A.copy()
    P = A.copy()
    for p in range(2, max(ks) + 1):
        P = ((P @ A) > 0).astype(np.int8)
        R = ((R + P) > 0).astype(np.int8)
        if p in ks:
            T = sp.triu(R, k=1).tocoo()
            out[p] = sorted(zip(T.row.tolist(), T.col.tolist()))
    return out


# ---------------------------------------------------------------- control
def null_edges(lat, lon, k):
    from scipy.spatial import cKDTree
    n = len(lat)
    X = xyz(lat, lon)
    tree = cKDTree(X)
    q = min(n, k + 12)
    _, idx = tree.query(X, k=q)
    E = set()
    for i in range(n):
        j = np.array([x for x in idx[i] if x != i])
        d = hav(lat[i], lon[i], lat[j], lon[j])
        o = np.lexsort((j, d))
        j, d = j[o], d[o]
        if q < n:   # the window must hold every node tied with the k-th
            dk = d[k - 1]
            assert d[-1] > dk + 1e-6, "neighbour window too small"
        for jj in j[:k]:
            jj = int(jj)
            E.add((min(i, jj), max(i, jj)))
    return sorted(E)


def round_half_up(x):
    return int(math.floor(x + 0.5))


# ---------------------------------------------------------------- stats
def graph_stats(n, edges, groups):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    deg = np.array([d for _, d in G.degree()])
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                degree_min=int(deg.min()), degree_median=float(np.median(deg)),
                degree_max=int(deg.max()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                frac_edges_within_group=float(np.mean(
                    [bool(memb[a] & memb[b]) for a, b in edges])))


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


def write_common(od, groups, lat, lon):
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
        for g in groups:
            fh.write(" ".join(map(str, g)) + "\n")
    with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
        fh.write("id\tlat\tlon\n")
        for i in range(len(lat)):
            fh.write(f"{i}\t{lat[i]:.6f}\t{lon[i]:.6f}\n")


# ---------------------------------------------------------------- build
def build_city(c, dl, out):
    path = os.path.join(dl, c["city"], "feed.zip")
    info = dict(city=c["city"], agency=c["agency"], url=c["url"],
                licence=c["licence"], licence_url=c["licence_url"],
                feed_bytes=os.path.getsize(path), feed_sha256=sha256(path),
                centre=list(c["centre"]))
    zf = zipfile.ZipFile(path)
    info["feed_info"] = feed_info(zf)
    routes = read(zf, "routes.txt", ["route_id", "route_type"])
    routes["rt"] = pd.to_numeric(routes["route_type"], errors="coerce").fillna(-1).astype(int)
    rail_r = set(routes.loc[routes["rt"].map(is_rail), "route_id"])
    bus_r = set(routes.loc[routes["rt"].map(is_bus), "route_id"])
    info["route_types_in_feed"] = {str(k): int(v) for k, v in
                                   sorted(Counter(routes["rt"]).items())}
    trips = read(zf, "trips.txt", ["route_id", "trip_id"])
    stops = read(zf, "stops.txt", ["stop_id", "stop_name", "stop_lat", "stop_lon",
                                   "location_type", "parent_station"])

    def load(rset):
        tr = trips[trips["route_id"].isin(rset)].reset_index(drop=True)
        code = pd.Series(np.arange(len(tr)), index=tr["trip_id"].to_numpy())
        assert code.index.is_unique
        st = read_stop_times(zf, code)
        return tr, st

    rtr, rst = load(rail_r)
    smap, pos, mstat = station_map(stops, set(rst["s"]))
    rpats = patterns(rst, rtr["route_id"].to_numpy(), smap)
    n_rail = len({x for _, seq in rpats if len(seq) >= 2 for x in seq})
    info.update(n_rail_routes=len(set(rtr["route_id"])), n_rail_trips=len(rtr),
                n_rail_stations=n_rail)
    mode = "rail" if n_rail >= MIN_RAIL_N else "railbus"
    pats = {(r, s, False) for r, s in rpats}
    if mode == "railbus":
        btr, bst = load(bus_r)
        st_all = pd.concat([rst, bst.assign(t=bst["t"] + len(rtr))], ignore_index=True)
        smap, pos, mstat = station_map(stops, set(st_all["s"]))
        rpats = patterns(rst, rtr["route_id"].to_numpy(), smap)
        bpats = patterns(bst, btr["route_id"].to_numpy(), smap)
        pats = {(r, s, False) for r, s in rpats} | {(r, s, True) for r, s in bpats}
        info.update(n_bus_routes=len(set(btr["route_id"])), n_bus_trips=len(btr))
        del st_all
    info["mode"] = mode
    info.update(mstat)
    inside = {}

    def clip(s):
        if s not in inside:
            la, lo = pos[s]
            inside[s] = bool(in_box(np.array([la]), np.array([lo]), c["centre"])[0])
        return inside[s]
    E1, route_st = l_space(pats, clip)
    keys = sorted({x for e in E1 for x in e})
    idx = {k: i for i, k in enumerate(keys)}
    n = len(keys)
    # rounded as written to nodes.tsv, so `null` rebuilds identical controls
    lat = np.round(np.array([pos[k][0] for k in keys]), 6)
    lon = np.round(np.array([pos[k][1] for k in keys]), 6)
    assert np.isfinite(lat).all() and np.isfinite(lon).all()
    E1i = sorted((min(idx[a], idx[b]), max(idx[a], idx[b])) for a, b in E1)
    rtype = dict(zip(routes["route_id"], routes["rt"]))
    # groups
    groups, seen, n_small, n_dup = [], {}, 0, 0
    g_route = []
    for r in sorted(route_st):
        g = tuple(sorted(idx[s] for s in route_st[r]))
        if len(g) < MIN_GROUP:
            n_small += 1
            continue
        if g in seen:
            n_dup += 1
            continue
        seen[g] = r
        groups.append(g)
        g_route.append(r)
    sizes = np.array([len(g) for g in groups])
    cnt = Counter(u for g in groups for u in g)
    gpn = np.array([cnt[u] for u in range(n)])
    s1 = graph_stats(n, E1i, groups)
    info.update(n=n, n_routes_with_steps=len(route_st), n_groups=len(groups),
                n_groups_rail=int(sum(is_rail(rtype[r]) for r in g_route)),
                n_groups_bus=int(sum(is_bus(rtype[r]) for r in g_route)),
                n_groups_dropped_small=n_small, n_groups_dropped_dup=n_dup,
                group_size_min=int(sizes.min()), group_size_median=float(np.median(sizes)),
                group_size_mean=float(sizes.mean()), group_size_max=int(sizes.max()),
                frac_nodes_in_group=float(np.mean(gpn >= 1)),
                groups_per_node_mean=float(gpn.mean()),
                frac_nodes_ge2_groups=float(np.mean(gpn >= 2)),
                groups_per_node_max=int(gpn.max()),
                g1_m=s1["m"], g1_mean_degree=s1["mean_degree"],
                g1_degree_max=s1["degree_max"], g1_clustering=s1["clustering"],
                g1_n_components=s1["n_components"],
                g1_frac_edges_within_group=s1["frac_edges_within_group"])
    if mode == "railbus":
        info["n_nodes_in_box"] = int(sum(clip(k) for k in keys))
    rows = []
    P = powers(n, E1i)
    for p, E in P.items():
        name = f"{c['city']}_{mode}_r{p}"
        s_e = graph_stats(n, E, groups)
        k = round_half_up(s_e["degree_median"])
        N = null_edges(lat, lon, k)
        s_n = graph_stats(n, N, groups)
        for suffix, EE in (("", E), ("_null", N)):
            od = os.path.join(out, name + suffix)
            write_common(od, groups, lat, lon)
            write_edges(os.path.join(od, "edges.tsv"), EE)
        row = dict(graph=name, control=name + "_null", power=p, null_k=k,
                   **{x: info[x] for x in info if x not in ("feed_info",)},
                   feed_info=info["feed_info"], **s_e,
                   **{"null_" + key: v for key, v in s_n.items()})
        row["n"] = n
        rows.append(row)
        print(json.dumps({x: row[x] for x in ("graph", "n", "m", "mean_degree",
              "clustering", "n_groups", "groups_per_node_mean",
              "frac_edges_within_group", "null_k", "null_frac_edges_within_group")}),
              flush=True)
    return rows


def cmd_build(args):
    os.makedirs(args.out, exist_ok=True)
    want = args.cities.split(",") if args.cities else [c["city"] for c in CITIES]
    sp = os.path.join(args.out, "summary.json")
    old = []
    if args.cities and os.path.isfile(sp):
        old = [r for r in json.load(open(sp)) if r["city"] not in want]
    summary = list(old)
    for c in CITIES:
        if c["city"] in want:
            print("==", c["city"], flush=True)
            summary += build_city(c, args.dl, args.out)
    order = {c["city"]: i for i, c in enumerate(CITIES)}
    summary.sort(key=lambda r: (order[r["city"]], r["power"]))
    with open(sp, "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def cmd_null(args):
    """Rewrite only the *_null folders from the real folders' files."""
    for name in sorted(os.listdir(args.out)):
        src = os.path.join(args.out, name)
        if name.endswith("_null") or not os.path.isfile(os.path.join(src, "edges.tsv")):
            continue
        rows = [l.split("\t") for l in open(os.path.join(src, "nodes.tsv")).read().splitlines()[1:]]
        lat = np.array([float(x[1]) for x in rows])
        lon = np.array([float(x[2]) for x in rows])
        deg = np.zeros(len(rows), int)
        for l in open(os.path.join(src, "edges.tsv")):
            a, b = map(int, l.split())
            deg[a] += 1
            deg[b] += 1
        k = round_half_up(float(np.median(deg)))
        dst = os.path.join(args.out, name + "_null")
        os.makedirs(dst, exist_ok=True)
        N = null_edges(lat, lon, k)
        write_edges(os.path.join(dst, "edges.tsv"), N)
        for f in ("groups.txt", "nodes.tsv"):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
        print(name + "_null", len(rows), k, len(N))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("build")
    a.add_argument("--dl", required=True, help="folder holding <city>/feed.zip")
    a.add_argument("--out", required=True)
    a.add_argument("--cities", default="")
    a = sub.add_parser("null")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"build": cmd_build, "null": cmd_null}[args.cmd](args)


if __name__ == "__main__":
    main()
