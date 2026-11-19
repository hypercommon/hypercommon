"""Build small real networks from commuting flows between areal units, with
administrative or functional regions as ground-truth groups, plus a pure
geometry control for each graph.

Families (all written to paper/analysis/data_geo/<name>/):
  us_county_state_k{8,16}      US counties (contiguous states + DC), groups = states
  us_county_cz_k{8,16}         same graph, groups = USDA ERS 2020 commuting zones
  ew_ltla_region_k{8,16}       England and Wales local authorities (LTLA), groups = regions
  ew_msoa_north_lad_k{8,16}    MSOAs of NE, NW, Yorkshire and the Humber, East and
                               West Midlands, groups = local authority districts
  ew_msoa_south_lad_k{8,16}    MSOAs of East of England, London, SE, SW and Wales,
                               groups = local authority districts
  <name>_null                  same nodes and groups; edges = k nearest units by
                               great-circle distance between the units' points

Source files (downloaded by hand into separate folders under --raw; never executed):
  acs_flows/table1.xlsx                         ACS 2016-2020 county-to-county commuting flows
  gazetteer2022_x/2022_Gaz_counties_national.txt Census 2022 Gazetteer, counties
  ers_cz/2020-commuting-zones.csv               USDA ERS 2020 commuting zones
  ons_odwp01ew_x/ODWP01EW_LTLA.csv, ODWP01EW_MSOA.csv  Census 2021 ODWP01EW (Nomis bulk zip)
  ons_lookup/page_*.json                        ONS MSOA21 -> LAD22 -> RGN22 lookup (ArcGIS REST pages)
  ons_msoa_pts/page_*.json                      ONS MSOA 2021 boundary attributes (LAT, LONG, BNG_E, BNG_N)
URLs and checksums are in data_geo/README.md.

Usage (from paper/, with an isolated interpreter that has numpy, pandas,
openpyxl and networkx):
  python -I analysis/geo_build.py --raw D:/tmp/geo_raw --out analysis/data_geo

Edge rule (fixed before any graph was looked at; same for every family):
  1. F[a,b] = workers living in unit a and working in unit b, for a != b, both
     units in the node set (self-flow, flows to units outside the node set, abroad,
     offshore or "no fixed place" are ignored).
  2. Symmetrise: S = F + F^T.
  3. Coverage normalisation (the Hi-C "VC_SQRT" rule): R[a] = sum_b S[a,b];
     score[a,b] = S[a,b] / sqrt(R[a] * R[b]).
  4. Each unit a keeps the k units b != a with the largest score (score > 0 only;
     ties broken by smaller great-circle distance, then smaller id). The
     undirected, unweighted edge set is the union over a.
Control (<name>_null): each unit keeps its k nearest units by great-circle
distance between the units' points (ties by smaller id), union over units.
Diagnostics only (not written as graphs): the same top-k rule on the raw
symmetrised flow S (no normalisation).
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KS = (8, 16)
EARTH_R_KM = 6371.0088
NON_CONTIGUOUS = {"AK", "HI", "PR"}
EW_NORTH = ("E12000001", "E12000002", "E12000003", "E12000004", "E12000005")
EW_SOUTH = ("E12000006", "E12000007", "E12000008", "E12000009", "W92000004")


# ----------------------------------------------------------------------------
# generic graph pieces
# ----------------------------------------------------------------------------
def haversine_matrix(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    dla = la[:, None] - la[None, :]
    dlo = lo[:, None] - lo[None, :]
    h = np.sin(dla / 2) ** 2 + np.cos(la)[:, None] * np.cos(la)[None, :] * np.sin(dlo / 2) ** 2
    return 2 * EARTH_R_KM * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def topk_union(score, dist, k):
    """Union over rows of each row's k highest-score columns (score > 0),
    ties by smaller distance then smaller index."""
    n = score.shape[0]
    edges = set()
    for i in range(n):
        s = score[i]
        cand = np.nonzero(s > 0)[0]
        cand = cand[cand != i]
        if len(cand) == 0:
            continue
        order = np.lexsort((cand, dist[i, cand], -s[cand]))
        for j in cand[order[:k]]:
            j = int(j)
            edges.add((i, j) if i < j else (j, i))
    return sorted(edges)


def nearest_union(dist, k):
    """Control: each node's k nearest nodes by distance, ties by smaller id."""
    n = dist.shape[0]
    idx = np.arange(n)
    edges = set()
    for i in range(n):
        cand = idx[idx != i]
        order = np.lexsort((cand, dist[i, cand]))
        for j in cand[order[:k]]:
            j = int(j)
            edges.add((i, j) if i < j else (j, i))
    return sorted(edges)


def flow_scores(F):
    """F: dense n x n directed flow with zero diagonal. Returns (S, VC score)."""
    S = F + F.T
    np.fill_diagonal(S, 0.0)
    R = S.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(R > 0, 1.0 / np.sqrt(R), 0.0)
    return S, S * inv[:, None] * inv[None, :]


def write_graph(od, edges, groups, nodes_df):
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "edges.tsv"), "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)
    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
        for g in groups:
            fh.write(" ".join(map(str, g)) + "\n")
    nodes_df.to_csv(os.path.join(od, "nodes.tsv"), sep="\t", index=False,
                    lineterminator="\n", encoding="utf-8")


def make_groups(labels, order_key=None):
    """labels: per-node group label (None = no group). Groups of size < 3 dropped.
    Returns (groups, n_dropped_small, n_labels)."""
    by = {}
    for i, g in enumerate(labels):
        if g is None or (isinstance(g, float) and np.isnan(g)):
            continue
        by.setdefault(g, []).append(i)
    keys = sorted(by, key=order_key)
    groups = [sorted(by[g]) for g in keys if len(by[g]) >= 3]
    return groups, sum(len(by[g]) < 3 for g in keys), len(keys)


def graph_stats(n, edges, groups, dist):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    deg = np.array([d for _, d in G.degree()])
    el = np.array([dist[a, b] for a, b in edges])
    within = float(np.mean([bool(memb[a] & memb[b]) for a, b in edges]))
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                max_degree=int(deg.max()), min_degree=int(deg.min()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                median_edge_km=float(np.median(el)),
                frac_edges_within_group=within)


def group_stats(n, groups):
    sizes = [len(g) for g in groups]
    covered = set().union(*map(set, groups)) if groups else set()
    return dict(n_groups=len(groups), group_size_min=min(sizes),
                group_size_median=float(np.median(sizes)),
                group_size_max=max(sizes), frac_nodes_in_group=len(covered) / n)


def build_family(out, base, nodes_df, F, lat, lon, groupsets, extra, summary):
    """groupsets: list of (suffix, groups, group_info). Writes graphs + controls."""
    n = len(nodes_df)
    dist = haversine_matrix(lat, lon)
    S, vc = flow_scores(F)
    for k in KS:
        edges = topk_union(vc, dist, k)
        raw_edges = topk_union(S, dist, k)
        null_edges = nearest_union(dist, k)
        for suffix, groups, ginfo in groupsets:
            name = f"{base}_{suffix}_k{k}"
            write_graph(os.path.join(out, name), edges, groups, nodes_df)
            write_graph(os.path.join(out, name + "_null"), null_edges, groups, nodes_df)
            row = dict(graph=name, control=name + "_null", k=k, n=n)
            row.update(graph_stats(n, edges, groups, dist))
            row.update(group_stats(n, groups))
            row.update(ginfo)
            ns = graph_stats(n, null_edges, groups, dist)
            row.update({f"null_{key}": v for key, v in ns.items()})
            rs = graph_stats(n, raw_edges, groups, dist)
            row.update({f"rawflow_rule_{key}": rs[key] for key in
                        ("m", "max_degree", "clustering", "median_edge_km",
                         "frac_edges_within_group")})
            row["frac_edges_shared_with_null"] = len(set(edges) & set(null_edges)) / len(edges)
            row["n_units_zero_flow"] = int((S.sum(1) == 0).sum())
            row.update(extra)
            summary.append(row)
            print(json.dumps({x: row[x] for x in ("graph", "n", "m", "mean_degree",
                                                   "clustering", "n_groups",
                                                   "frac_edges_within_group",
                                                   "null_frac_edges_within_group",
                                                   "n_components")}))


# ----------------------------------------------------------------------------
# United States
# ----------------------------------------------------------------------------
def build_us(raw, out, summary):
    gaz = pd.read_csv(os.path.join(raw, "gazetteer2022_x", "2022_Gaz_counties_national.txt"),
                      sep="\t", dtype=str)
    gaz.columns = [c.strip() for c in gaz.columns]
    gaz = gaz[~gaz["USPS"].isin(NON_CONTIGUOUS)].copy()
    gaz = gaz.sort_values("GEOID").reset_index(drop=True)
    fips = gaz["GEOID"].tolist()
    pos = {f: i for i, f in enumerate(fips)}
    n = len(fips)

    cz = pd.read_csv(os.path.join(raw, "ers_cz", "2020-commuting-zones.csv"), dtype=str)
    cz = cz.set_index("FIPStxt")
    missing_cz = [f for f in fips if f not in cz.index]
    assert not missing_cz, missing_cz[:10]

    fl = pd.read_excel(os.path.join(raw, "acs_flows", "table1.xlsx"), header=None,
                       skiprows=8, dtype=str)
    fl = fl[fl[0].str.fullmatch(r"\d{2}", na=False) & fl[1].str.fullmatch(r"\d{3}", na=False)]
    # workplace state "0SS" = US state SS; 301/303 (Canada/Mexico) and the
    # "outside of the U.S." rows have no county and are dropped
    fl = fl[fl[4].str.fullmatch(r"0\d{2}", na=False) & fl[5].str.fullmatch(r"\d{3}", na=False)]
    orig = fl[0] + fl[1]
    dest = fl[4].str[1:] + fl[5]
    cnt = fl[8].astype(int).to_numpy()
    res_counties = set(orig)
    missing_acs = [f for f in fips if f not in res_counties]
    assert not missing_acs, missing_acs[:10]
    F = np.zeros((n, n))
    oi = orig.map(pos)
    di = dest.map(pos)
    ok = oi.notna() & di.notna()
    oi, di, c = oi[ok].astype(int).to_numpy(), di[ok].astype(int).to_numpy(), cnt[ok.to_numpy()]
    off = oi != di
    np.add.at(F, (oi[off], di[off]), c[off])
    extra = dict(country="US", unit="county (2022 Gazetteer; CT planning regions)",
                 flow_source="ACS 2016-2020 5-year commuting flows, table 1",
                 flow_rows_used=int(off.sum()), flow_workers_cross=float(F.sum()),
                 flow_workers_self=float(c[~off].sum()))

    nodes = pd.DataFrame(dict(
        id=np.arange(n), fips=fips, name=gaz["NAME"].str.strip(), state=gaz["USPS"],
        lat=gaz["INTPTLAT"].astype(float),
        lon=gaz["INTPTLONG"].str.strip().astype(float),
        cz2020=[cz.loc[f, "CZ2020"] for f in fips],
        cz2020_name=[cz.loc[f, "CZName"] for f in fips]))
    st_groups, st_small, st_all = make_groups([f[:2] for f in fips])
    cz_groups, cz_small, cz_all = make_groups(nodes["cz2020"].tolist(), order_key=int)
    groupsets = [
        ("state", st_groups, dict(group_source="state (+DC)", n_labels=st_all,
                                  n_groups_dropped_small=st_small)),
        ("cz", cz_groups, dict(group_source="USDA ERS 2020 commuting zones",
                               n_labels=cz_all, n_groups_dropped_small=cz_small)),
    ]
    build_family(out, "us_county", nodes, F, nodes["lat"].to_numpy(),
                 nodes["lon"].to_numpy(), groupsets, extra, summary)


# ----------------------------------------------------------------------------
# England and Wales
# ----------------------------------------------------------------------------
def read_pages(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)):
        d = json.load(open(f, encoding="utf-8"))
        assert "error" not in d, d
        rows += [x["attributes"] for x in d["features"]]
    return pd.DataFrame(rows)


def od_matrix(path, code_col, dest_col, pos):
    od = pd.read_csv(path, usecols=[code_col, dest_col, "Count"], dtype={code_col: str, dest_col: str})
    oi = od[code_col].map(pos)
    di = od[dest_col].map(pos)
    ok = oi.notna() & di.notna()
    oi, di = oi[ok].astype(int).to_numpy(), di[ok].astype(int).to_numpy()
    c = od.loc[ok, "Count"].to_numpy().astype(float)
    n = len(pos)
    F = np.zeros((n, n))
    off = oi != di
    np.add.at(F, (oi[off], di[off]), c[off])
    return F, int(off.sum()), float(c[~off].sum())


def build_ew(raw, out, summary):
    lu = read_pages(os.path.join(raw, "ons_lookup", "page_*.json"))
    pts = read_pages(os.path.join(raw, "ons_msoa_pts", "page_*.json"))
    assert lu["MSOA21CD"].is_unique and pts["MSOA21CD"].is_unique and len(lu) == len(pts) == 7264
    ms = lu.merge(pts[["MSOA21CD", "BNG_E", "BNG_N", "LAT", "LONG"]], on="MSOA21CD",
                  validate="1:1").sort_values("MSOA21CD").reset_index(drop=True)
    odir = os.path.join(raw, "ons_odwp01ew_x")
    msoa_csv = os.path.join(odir, "ODWP01EW_MSOA.csv")
    ltla_csv = os.path.join(odir, "ODWP01EW_LTLA.csv")
    src = "Census 2021 ODWP01EW (Nomis bulk)"

    # resident workers per MSOA (all rows, any workplace) - only for LTLA points
    od = pd.read_csv(msoa_csv, usecols=["Middle layer Super Output Areas code", "Count"],
                     dtype=str)
    workers = od.groupby("Middle layer Super Output Areas code")["Count"].apply(
        lambda s: s.astype(int).sum())
    ms["workers"] = ms["MSOA21CD"].map(workers).astype(float)
    assert ms["workers"].notna().all()

    # --- LTLA level, groups = regions
    lad = (ms.groupby(["LAD22CD", "LAD22NM", "RGN22CD", "RGN22NM"], as_index=False)
           .apply(lambda d: pd.Series(dict(
               lat=float(np.average(d["LAT"], weights=d["workers"])),
               lon=float(np.average(d["LONG"], weights=d["workers"])),
               n_msoa=len(d))), include_groups=False)
           .sort_values("LAD22CD").reset_index(drop=True))
    assert lad["LAD22CD"].is_unique and len(lad) == 331
    pos = {c: i for i, c in enumerate(lad["LAD22CD"])}
    F, nrows, selfw = od_matrix(ltla_csv, "Lower tier local authorities code",
                                "LTLA of workplace code", pos)
    nodes = pd.DataFrame(dict(id=np.arange(len(lad)), lad22cd=lad["LAD22CD"],
                              name=lad["LAD22NM"], rgn22cd=lad["RGN22CD"],
                              rgn22nm=lad["RGN22NM"], lat=lad["lat"].round(5),
                              lon=lad["lon"].round(5), n_msoa=lad["n_msoa"].astype(int)))
    groups, small, nl = make_groups(nodes["rgn22cd"].tolist())
    extra = dict(country="England and Wales", unit="LTLA / LAD 2022",
                 flow_source=src + ", LTLA table", flow_rows_used=nrows,
                 flow_workers_cross=float(F.sum()), flow_workers_self=selfw,
                 point="mean of the LTLA's MSOA points weighted by resident workers")
    build_family(out, "ew_ltla", nodes, F, nodes["lat"].to_numpy(), nodes["lon"].to_numpy(),
                 [("region", groups, dict(group_source="ONS region (RGN22; Wales as one)",
                                          n_labels=nl, n_groups_dropped_small=small))],
                 extra, summary)

    # --- MSOA level, two halves, groups = LADs
    for part, rgns in (("north", EW_NORTH), ("south", EW_SOUTH)):
        sub = ms[ms["RGN22CD"].isin(rgns)].reset_index(drop=True)
        pos = {c: i for i, c in enumerate(sub["MSOA21CD"])}
        F, nrows, selfw = od_matrix(msoa_csv, "Middle layer Super Output Areas code",
                                    "MSOA of workplace code", pos)
        nodes = pd.DataFrame(dict(id=np.arange(len(sub)), msoa21cd=sub["MSOA21CD"],
                                  name=sub["MSOA21NM"], lad22cd=sub["LAD22CD"],
                                  lad22nm=sub["LAD22NM"], rgn22nm=sub["RGN22NM"],
                                  lat=sub["LAT"], lon=sub["LONG"],
                                  bng_e=sub["BNG_E"], bng_n=sub["BNG_N"]))
        groups, small, nl = make_groups(nodes["lad22cd"].tolist())
        extra = dict(country="England and Wales", unit=f"MSOA 2021 ({part}: {', '.join(rgns)})",
                     flow_source=src + ", MSOA table (flows within the part only)",
                     flow_rows_used=nrows, flow_workers_cross=float(F.sum()),
                     flow_workers_self=selfw,
                     point="ONS MSOA 2021 boundary LAT/LONG (BSC V3 attributes)")
        build_family(out, f"ew_msoa_{part}", nodes, F, nodes["lat"].to_numpy(),
                     nodes["lon"].to_numpy(),
                     [("lad", groups, dict(group_source="local authority district (LAD22)",
                                           n_labels=nl, n_groups_dropped_small=small))],
                     extra, summary)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, help="folder holding the download folders")
    ap.add_argument("--out", required=True, help="paper/analysis/data_geo")
    ap.add_argument("--only", choices=["us", "ew"], default=None)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    summary = []
    if args.only in (None, "us"):
        build_us(args.raw, args.out, summary)
    if args.only in (None, "ew"):
        build_ew(args.raw, args.out, summary)
    if args.only is None:
        with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
            json.dump(summary, fh, indent=1)


if __name__ == "__main__":
    main()
