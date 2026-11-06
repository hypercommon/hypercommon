"""
Loaders for the additional real networks with recorded groups, from the
Netzschleuder repository (https://networks.skewed.de), files in data_netz/.

Each loader returns (G, truth): G a simple undirected networkx graph (self-loops
and directions dropped), truth a list of node sets with at least 3 members.

  facebook_<ego>   Facebook ego networks (McAuley & Leskovec 2012): the friends of
                   one user; ground truth = the user's own social circles, which
                   overlap. The ego itself is not in the files.
  cora             citation network; ground truth = paper topic (partition)
  polblogs         political blogs, 2004; ground truth = liberal / conservative
  primary_school   face-to-face contacts in a primary school (SocioPatterns);
                   ground truth = school class (teachers excluded)
  football         US college football 2000; ground truth = conference
  polbooks         co-purchased political books; ground truth = leaning
"""
from __future__ import annotations

import ast
import csv
import io
import os
import re
import zipfile

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
NETZ = os.path.join(HERE, "data_netz")
FACEBOOK_EGOS = [0, 107, 348, 414, 686, 698, 1684, 1912, 3437, 3980]


def _read(zipname):
    with zipfile.ZipFile(os.path.join(NETZ, zipname)) as z:
        nodes = z.read("nodes.csv").decode("utf-8")
        edges = z.read("edges.csv").decode("utf-8")
    nrows = list(csv.reader(io.StringIO(nodes)))
    header = [h.strip().lstrip("#").strip() for h in nrows[0]]
    nodes = [dict(zip(header, [c.strip() for c in r])) for r in nrows[1:] if r]
    erows = list(csv.reader(io.StringIO(edges)))
    G = nx.Graph()
    G.add_nodes_from(int(r["index"]) for r in nodes)
    for r in erows[1:]:
        if r and not r[0].startswith("#"):
            a, b = int(r[0]), int(r[1])
            if a != b:
                G.add_edge(a, b)
    return G, nodes


def _array(s):
    """Parse numpy-style text such as 'array([1, 5], dtype=int32)' or '[1 5]'."""
    nums = re.findall(r"-?\d+", s.split("dtype")[0])
    return [int(x) for x in nums]


def _groups_from_label(nodes, key, skip=()):
    groups = {}
    for r in nodes:
        v = r.get(key, "")
        if v == "" or v in skip:
            continue
        groups.setdefault(v, set()).add(int(r["index"]))
    return [g for g in groups.values() if len(g) >= 3]


def load(name):
    if name.startswith("facebook_"):
        G, nodes = _read(f"ego_social__{name}.csv.zip")
        circles = {}
        for r in nodes:
            for c in _array(r.get("circles", "")):
                circles.setdefault(c, set()).add(int(r["index"]))
        truth = [c for c in circles.values() if len(c) >= 3]
        return G, truth
    if name == "cora":
        G, nodes = _read("cora__cora.csv.zip")
        return G, _groups_from_label(nodes, "meta")
    if name == "polblogs":
        G, nodes = _read("polblogs__polblogs.csv.zip")
        return G, _groups_from_label(nodes, "value")
    if name == "primary_school":
        G, nodes = _read("sp_primary_school__sp_primary_school.csv.zip")
        return G, _groups_from_label(nodes, "class", skip=("Teachers",))
    if name == "football":
        G, nodes = _read("football__football.csv.zip")
        return G, _groups_from_label(nodes, "value")
    if name == "polbooks":
        G, nodes = _read("polbooks__polbooks.csv.zip")
        return G, _groups_from_label(nodes, "value")
    raise ValueError(name)


ALL = [f"facebook_{e}" for e in FACEBOOK_EGOS] + ["cora", "polblogs", "primary_school", "football", "polbooks"]


if __name__ == "__main__":
    for nm in ALL:
        G, truth = load(nm)
        covered = set().union(*truth) if truth else set()
        memb = sum(len(c) for c in truth)
        print(f"{nm:16s} n={G.number_of_nodes():6d} m={G.number_of_edges():7d} clustering={nx.transitivity(G):.3f} "
              f"groups={len(truth):3d} labelled={len(covered)/G.number_of_nodes():.2f} "
              f"memberships/labelled={memb/max(len(covered),1):.2f}")
