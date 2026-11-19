"""Build street graphs of city centres from OpenStreetMap, with named streets as
ground-truth groups. Groups are long thin chains that cross and share nodes:
a point at a crossing belongs to every street through it.

Source data (fetched by the `fetch` command, one Overpass request per city,
saved untouched into its own folder and only ever parsed as JSON):
  * OpenStreetMap, (c) OpenStreetMap contributors, ODbL 1.0,
    via the public Overpass API https://overpass-api.de/api/interpreter
    (query text, sha256 and OSM timestamps in data_streets/README.md).

Usage (the inputs are untrusted downloads, so run with an isolated interpreter):
  python -I analysis/streets_build.py fetch --dl <empty dir>
  python -I analysis/streets_build.py build --dl <same dir> --out analysis/data_streets
(run from paper/; any working directory works, all paths are arguments).

Rules (fixed before any graph was built or any statistic looked at; no
community detection is run):
  Areas: the 12 fixed centres in CITIES. Each city has one download: a square
    box of side OUTER_SIDE_M around the centre (big enough for every box below).
    Boxes are squares in a local equirectangular projection around the centre
    (x = R cos(lat0) dlon, y = R dlat, R = 6371008.8 m).
  Box size rule: sides s_j = 1500 m * 1.25**j, j in -4..4. Start at j = 0;
    while n (below) > 5000 take j - 1, while n < 2000 take j + 1; stop at the
    first j with 2000 <= n <= 5000. If j would leave -4..4 the city is
    dropped (recorded). If the walk turns round (n jumps over the range) the
    j with n closest to the range in log(n) is used (recorded).
  Ways: highway in HIGHWAYS with a `name` tag; ways tagged area=yes
    (pedestrian squares drawn as outlines) are excluded.
  Clip: a way segment (two consecutive OSM nodes) is kept iff both its end
    nodes lie inside the box. Segments joining the same node pair are merged.
  Nodes: OSM nodes that end a kept segment, plus interior points on each kept
    segment: a segment of length L gets ceil(L / 30 m) - 1 evenly spaced
    interior points (so consecutive points are <= 30 m apart). Interior
    points are keyed by (smaller OSM id, larger OSM id, index from the smaller
    end), so a segment shared by two ways gets one set of points. A crossing
    (a shared OSM node) is one node. Positions are linear in lat/lon.
    Node ids: OSM nodes in OSM-id order, then interior points in key order.
  G1: consecutive points along each kept segment are linked.
  Graphs: <city>_r2 = G1 squared (pairs within 2 hops of G1),
          <city>_r3 = G1 cubed (within 3 hops). No self-loops.
  Control <name>_null: same nodes and groups; k = floor(median degree of the
    real graph); each node keeps its k nearest other nodes by straight-line
    (projected) distance, ties by smaller id; union over nodes.
  Groups: all kept segments of ways with the same exact `name` string are
    merged; the group is the set of their points (end nodes and interior
    points). Groups with < 3 points are dropped; an exact duplicate point set
    is kept once (first in name order). Groups are written in code-point
    order of the name (names themselves are not written).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

R_EARTH = 6371008.8
OUTER_SIDE_M = 3700.0
BASE_SIDE_M = 1500.0
SIDE_FACTOR = 1.25
J_RANGE = (-4, 4)
N_MIN, N_MAX = 2000, 5000
SPACING_M = 30.0
MIN_GROUP = 3
HIGHWAYS = ("primary", "secondary", "tertiary", "residential", "unclassified",
            "living_street", "pedestrian")
OSM_DATE = "2026-11-15T00:00:00Z"
ENDPOINT = "https://overpass-api.de/api/interpreter"
USER_AGENT = "hypercommon-benchmark/1.0 (research graph build; python-urllib)"

# name: (lat, lon, description). Fixed before any download.
CITIES = {
    "manhattan_midtown": (40.7549, -73.9840, "Manhattan, Bryant Park / Times Square"),
    "chicago_loop": (41.8807, -87.6298, "Chicago, the Loop"),
    "barcelona_eixample": (41.3917, 2.1649, "Barcelona, Eixample"),
    "paris_centre": (48.8589, 2.3469, "Paris, Chatelet"),
    "rome_centre": (41.8986, 12.4769, "Rome, Pantheon"),
    "tokyo_ginza": (35.6717, 139.7650, "Tokyo, Ginza"),
    "kyoto_centre": (35.0060, 135.7620, "Kyoto, Shijo-Karasuma"),
    "amsterdam_centre": (52.3731, 4.8926, "Amsterdam, Dam"),
    "london_centre": (51.5130, -0.1300, "London, Soho / Covent Garden"),
    "prague_old_town": (50.0870, 14.4210, "Prague, Old Town"),
    "melbourne_cbd": (-37.8136, 144.9631, "Melbourne, CBD"),
    "buenos_aires_centro": (-34.6037, -58.3816, "Buenos Aires, Obelisco"),
}


# ---------------------------------------------------------------- geometry
def box_deg(lat0, lon0, side_m):
    h = side_m / 2.0
    dlat = math.degrees(h / R_EARTH)
    dlon = math.degrees(h / (R_EARTH * math.cos(math.radians(lat0))))
    return lat0 - dlat, lon0 - dlon, lat0 + dlat, lon0 + dlon   # s, w, n, e


def project(lat, lon, lat0, lon0):
    x = R_EARTH * math.cos(math.radians(lat0)) * math.radians(lon - lon0)
    y = R_EARTH * math.radians(lat - lat0)
    return x, y


# ---------------------------------------------------------------- fetch
def overpass_query(city):
    lat0, lon0, _ = CITIES[city]
    s, w, n, e = box_deg(lat0, lon0, OUTER_SIDE_M)
    hw = "|".join(HIGHWAYS)
    return (f'[out:json][timeout:300][date:"{OSM_DATE}"];\n'
            f'way["highway"~"^({hw})$"]["name"]'
            f"({s:.7f},{w:.7f},{n:.7f},{e:.7f});\n"
            f"out body geom;\n")


def post(query, tries=6):
    data = ("data=" + urllib.parse.quote(query)).encode()
    wait = 30
    for t in range(tries):
        req = urllib.request.Request(ENDPOINT, data=data, headers={
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded"})
        try:
            with urllib.request.urlopen(req, timeout=400) as r:
                body = r.read()
            js = json.loads(body)
            remark = js.get("remark", "")
            if "error" in remark.lower() or "timed out" in remark.lower():
                raise RuntimeError("overpass remark: " + remark)
            return body
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError,
                RuntimeError, json.JSONDecodeError, ConnectionError) as ex:
            print(f"  try {t + 1} failed: {ex!r}; waiting {wait} s", flush=True)
            time.sleep(wait)
            wait *= 2
    raise SystemExit("giving up")


def cmd_fetch(args):
    os.makedirs(args.dl, exist_ok=True)
    for city in CITIES:
        jf = os.path.join(args.dl, city + ".json")
        qf = os.path.join(args.dl, city + ".overpassql")
        q = overpass_query(city)
        if os.path.isfile(jf) and os.path.isfile(qf) and open(qf).read() == q:
            print(city, "already there")
            continue
        print(city, "fetching", flush=True)
        body = post(q)
        with open(qf, "w", newline="\n") as fh:
            fh.write(q)
        with open(jf, "wb") as fh:
            fh.write(body)
        print(city, len(body), hashlib.sha256(body).hexdigest(), flush=True)
        time.sleep(15)                       # one request at a time, gently


# ---------------------------------------------------------------- build
def load_ways(dl, city):
    with open(os.path.join(dl, city + ".json"), "rb") as fh:
        raw = fh.read()
    js = json.loads(raw)
    ways = []
    for el in js["elements"]:
        if el.get("type") != "way":
            continue
        tags = el.get("tags", {})
        if tags.get("highway") not in HIGHWAYS or "name" not in tags:
            continue
        if tags.get("area") == "yes":
            continue
        geom = el.get("geometry")
        nodes = el.get("nodes")
        assert geom is not None and len(geom) == len(nodes), el["id"]
        ways.append(dict(id=el["id"], name=tags["name"], nodes=nodes,
                         coords=[(g["lat"], g["lon"]) for g in geom]))
    meta = dict(osm3s=js.get("osm3s", {}), generator=js.get("generator"),
                sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw),
                n_ways_raw=sum(1 for el in js["elements"] if el.get("type") == "way"))
    return ways, meta


def clip_points(ways, lat0, lon0, side_m):
    """Return node coords, G1 edges, group dict name -> set of node ids, info."""
    h = side_m / 2.0
    pos = {}                                    # osm id -> (lat, lon)
    segs = {}                                   # (a, b) a < b -> set of names
    for w in ways:
        for (a, b), (ca, cb) in zip(zip(w["nodes"], w["nodes"][1:]),
                                    zip(w["coords"], w["coords"][1:])):
            if a == b:
                continue
            xa, ya = project(*ca, lat0, lon0)
            xb, yb = project(*cb, lat0, lon0)
            if max(abs(xa), abs(ya), abs(xb), abs(yb)) > h:
                continue
            pos[a] = ca
            pos[b] = cb
            key = (a, b) if a < b else (b, a)
            segs.setdefault(key, set()).add(w["name"])
    osm_ids = sorted(pos)
    idx = {o: i for i, o in enumerate(osm_ids)}
    coords = [pos[o] for o in osm_ids]
    interior = []                               # (a, b, i, lat, lon)
    chain = {}                                  # seg -> list of node ids a..b
    total_len = 0.0
    for (a, b) in sorted(segs):
        (la, oa), (lb, ob) = pos[a], pos[b]
        xa, ya = project(la, oa, lat0, lon0)
        xb, yb = project(lb, ob, lat0, lon0)
        L = math.hypot(xb - xa, yb - ya)
        total_len += L
        q = max(1, math.ceil(L / SPACING_M))
        chain[(a, b)] = [(a, b, i, la + (lb - la) * i / q, oa + (ob - oa) * i / q)
                         for i in range(1, q)]
        interior.extend(chain[(a, b)])
    base = len(coords)
    iid = {}
    for k, p in enumerate(interior):            # already in key order
        iid[p[:3]] = base + k
        coords.append((p[3], p[4]))
    edges = set()
    groups = defaultdict(set)
    for (a, b), names in segs.items():
        path = [idx[a]] + [iid[p[:3]] for p in chain[(a, b)]] + [idx[b]]
        for u, v in zip(path, path[1:]):
            edges.add((min(u, v), max(u, v)))
        for nm in names:
            groups[nm].update(path)
    info = dict(n_osm_nodes=base, n_interior_points=len(interior),
                n_segments=len(segs), street_length_km=total_len / 1000.0,
                n_segments_multi_name=sum(1 for s in segs.values() if len(s) > 1))
    return coords, sorted(edges), groups, info


def choose_side(ways, lat0, lon0):
    j, seen = 0, {}
    while True:
        side = BASE_SIDE_M * SIDE_FACTOR ** j
        n = len(clip_points(ways, lat0, lon0, side)[0])
        seen[j] = n
        if N_MIN <= n <= N_MAX:
            return j, side, seen, "in range"
        nj = j - 1 if n > N_MAX else j + 1
        if nj in seen:                          # turned round
            best = min(seen, key=lambda t: min(abs(math.log(seen[t] / N_MIN)),
                                               abs(math.log(seen[t] / N_MAX))))
            return best, BASE_SIDE_M * SIDE_FACTOR ** best, seen, "turned round"
        if not (J_RANGE[0] <= nj <= J_RANGE[1]):
            return None, None, seen, "out of j range"
        j = nj


def power_edges(n, g1, r):
    import numpy as np
    import scipy.sparse as sp
    e = np.array(g1)
    A = sp.coo_matrix((np.ones(len(e) * 2, np.int64),
                       (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                      shape=(n, n)).tocsr()
    B = A + sp.identity(n, dtype=np.int64, format="csr")
    P = B
    for _ in range(r - 1):
        P = (P @ B).tocsr()
        P.data[:] = 1
    P = sp.triu(P, k=1).tocoo()
    return sorted(zip(P.row.tolist(), P.col.tolist()))


def null_edges(xy, k):
    import numpy as np
    from scipy.spatial import cKDTree
    n = len(xy)
    tree = cKDTree(xy)
    q = min(n, k + 17)
    d, j = tree.query(xy, k=q)
    E = set()
    for i in range(n):
        di, ji = d[i], j[i]
        keep = ji != i
        di, ji = di[keep], ji[keep]
        order = np.lexsort((ji, di))
        di, ji = di[order], ji[order]
        assert len(di) > k and di[k - 1] < d[i][-1], "query window too small"
        for v in ji[:k]:
            v = int(v)
            E.add((min(i, v), max(i, v)))
    return sorted(E)


def graph_stats(n, edges, memb, xy):
    import networkx as nx
    import numpy as np
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    e = np.array(edges)
    deg = np.bincount(e.ravel(), minlength=n)
    within = float(np.mean([bool(memb[a] & memb[b]) for a, b in edges]))
    lens = np.hypot(*(xy[e[:, 0]] - xy[e[:, 1]]).T)
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                degree_min=int(deg.min()), degree_median=float(np.median(deg)),
                degree_max=int(deg.max()),
                clustering=nx.average_clustering(G),
                n_components=len(comps), largest_component=comps[0],
                edge_length_median_m=float(np.median(lens)),
                edge_length_max_m=float(lens.max()),
                frac_edges_within_group=within)


def write_edges(path, edges):
    with open(path, "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)


def cmd_build(args):
    import numpy as np
    os.makedirs(args.out, exist_ok=True)
    summary, cities = [], []
    for city, (lat0, lon0, desc) in CITIES.items():
        ways, meta = load_ways(args.dl, city)
        j, side, seen, how = choose_side(ways, lat0, lon0)
        crow = dict(city=city, description=desc, centre_lat=lat0, centre_lon=lon0,
                    n_ways_named=len(ways), size_walk={str(t): v for t, v in seen.items()},
                    size_rule=how, j=j, side_m=side, download=meta)
        cities.append(crow)
        print(city, how, seen, flush=True)
        if j is None:
            continue
        coords, g1, groups_d, info = clip_points(ways, lat0, lon0, side)
        n = len(coords)
        groups, seen_sets, n_small, n_dup = [], set(), 0, 0
        for nm in sorted(groups_d):
            mem = tuple(sorted(groups_d[nm]))
            if len(mem) < MIN_GROUP:
                n_small += 1
                continue
            if mem in seen_sets:
                n_dup += 1
                continue
            seen_sets.add(mem)
            groups.append(mem)
        memb = [set() for _ in range(n)]
        for gi, g in enumerate(groups):
            for u in g:
                memb[u].add(gi)
        cnt = np.array([len(s) for s in memb])
        sizes = np.array([len(g) for g in groups])
        xy = np.array([project(la, lo, lat0, lon0) for la, lo in coords])
        s, w, nn, e = box_deg(lat0, lon0, side)
        g1_stats = graph_stats(n, g1, memb, xy)
        base = dict(city=city, side_m=side, box_swne=[s, w, nn, e], n=n,
                    **info, g1_m=g1_stats["m"], g1_n_components=g1_stats["n_components"],
                    g1_largest_component=g1_stats["largest_component"],
                    g1_frac_edges_within_group=g1_stats["frac_edges_within_group"],
                    n_groups=len(groups), n_groups_dropped_small=n_small,
                    n_groups_dropped_dup=n_dup,
                    group_size_min=int(sizes.min()),
                    group_size_median=float(np.median(sizes)),
                    group_size_mean=float(sizes.mean()),
                    group_size_max=int(sizes.max()),
                    frac_nodes_in_group=float(np.mean(cnt >= 1)),
                    groups_per_node_mean=float(cnt.mean()),
                    frac_nodes_in_ge2_groups=float(np.mean(cnt >= 2)),
                    groups_per_node_max=int(cnt.max()))
        for r in (2, 3):
            name = f"{city}_r{r}"
            E = power_edges(n, g1, r)
            s_e = graph_stats(n, E, memb, xy)
            k = int(math.floor(s_e["degree_median"]))
            N = null_edges(xy, k)
            s_n = graph_stats(n, N, memb, xy)
            for suffix, edges in (("", E), ("_null", N)):
                od = os.path.join(args.out, name + suffix)
                os.makedirs(od, exist_ok=True)
                write_edges(os.path.join(od, "edges.tsv"), edges)
                with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
                    for g in groups:
                        fh.write(" ".join(map(str, g)) + "\n")
                with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
                    fh.write("id\tlat\tlon\n")
                    for i, (la, lo) in enumerate(coords):
                        fh.write(f"{i}\t{la:.7f}\t{lo:.7f}\n")
            row = dict(graph=name, control=name + "_null", r=r, null_k=k, **base,
                       size_rule=how, size_walk=crow["size_walk"],
                       download_sha256=meta["sha256"],
                       osm_base=meta["osm3s"].get("timestamp_osm_base"), **s_e,
                       **{"null_" + key: v for key, v in s_n.items()})
            summary.append(row)
            print(json.dumps({x: row[x] for x in (
                "graph", "n", "m", "mean_degree", "clustering", "n_components",
                "n_groups", "groups_per_node_mean", "frac_nodes_in_ge2_groups",
                "frac_edges_within_group", "null_k", "null_m",
                "null_frac_edges_within_group")}, ensure_ascii=False), flush=True)
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)
    with open(os.path.join(args.out, "cities.json"), "w", newline="\n") as fh:
        json.dump(cities, fh, indent=1, ensure_ascii=False)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("fetch")
    a.add_argument("--dl", required=True, help="own empty folder for the raw responses")
    a = sub.add_parser("build")
    a.add_argument("--dl", required=True)
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"fetch": cmd_fetch, "build": cmd_build}[args.cmd](args)


if __name__ == "__main__":
    main()
