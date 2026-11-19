"""Build residue contact networks from protein 3D structures, with CATH
domains as ground-truth groups.

Usage (downloads are untrusted data; run with -I):
    python -I protein_build.py RAW_DIR OUT_DIR

RAW_DIR must already hold the CATH v4.4.0 files
    cath-domain-boundaries-v4_4_0.txt
    cath-domain-list-v4_4_0.txt
from http://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/
mmCIF files are fetched one by one from https://files.rcsb.org/download/<ID>.cif.gz
into RAW_DIR/cif/ (cached; reruns do not download again).

Writes OUT_DIR/<chain>/ and OUT_DIR/<chain>_7A/ (edges.tsv, groups.txt,
nodes.tsv), OUT_DIR/selection_log.tsv and OUT_DIR/summary.json.
"""
import json
import os
import random
import sys
import time
import urllib.request
from collections import defaultdict

import gemmi
import networkx as nx
import numpy as np

SEED = 2026
N_TARGET = 60
MIN_LEN, MAX_LEN = 150, 800          # residues with a C-alpha atom
MAX_RES = 2.5                        # Angstrom, X-ray only
MAX_INTERNAL_MISSING = 5             # SEQRES residues without CA between first and last observed CA
MAX_TOTAL_MISSING_FRAC = 0.10        # all SEQRES residues without CA (incl. termini), as fraction of SEQRES
MIN_FRAC_IN_DOMAIN = 0.80            # CA residues inside a kept CATH domain, as fraction of n
MAX_TRIES_PER_COMBO = 3
MAX_DOWNLOADS = 600
CUTOFFS = (8.0, 7.0)
MIN_DOMAIN = 3
RCSB = "https://files.rcsb.org/download/{}.cif.gz"


# ----------------------------------------------------------------- CATH parsing
def parse_boundaries(path):
    """CathDomall format 2.0 -> {chain: [ [ (ch, s, si, ch, e, ei), ... ] per domain ]}"""
    out = {}
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            chain, nd, nf = t[0], int(t[1][1:]), int(t[2][1:])
            k = 3
            doms = []
            for _ in range(nd):
                nseg = int(t[k]); k += 1
                segs = []
                for _ in range(nseg):
                    c1, s, si, c2, e, ei = t[k:k + 6]; k += 6
                    segs.append((c1, int(s), "" if si == "-" else si,
                                 c2, int(e), "" if ei == "-" else ei))
                doms.append(segs)
            out[chain] = (doms, nf)
    return out


def parse_domain_list(path):
    """-> {domain_id: (superfamily 'C.A.T.H', length, resolution)}"""
    out = {}
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            out[t[0]] = (".".join(t[1:5]), int(t[10]), float(t[11]))
    return out


# ------------------------------------------------------------------- download
_last = [0.0]


def fetch(pdb, cifdir):
    path = os.path.join(cifdir, pdb.lower() + ".cif.gz")
    if os.path.exists(path):
        return path, False
    wait = 1.0 - (time.time() - _last[0])
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(RCSB.format(pdb.upper()),
                                 headers={"User-Agent": "research-script/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
    except Exception as exc:  # noqa: BLE001
        _last[0] = time.time()
        return None, str(exc)
    _last[0] = time.time()
    with open(path, "wb") as fh:
        fh.write(data)
    return path, True


# ------------------------------------------------------------------ structure
def load_chain(path, chain_id):
    doc = gemmi.cif.read(path)
    block = doc.sole_block()
    method = gemmi.cif.as_string(block.find_value("_exptl.method") or "")
    nmeth = len(block.find_values("_exptl.method"))
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    info = {"method": method, "n_methods": nmeth, "resolution": st.resolution}
    model = st[0]                                    # first model
    ch = model.find_chain(chain_id)                  # author chain id
    if ch is None:
        return info, None
    poly = ch.get_polymer()
    if len(poly) == 0:
        return info, None
    ent = st.get_entity_of(poly)
    info["n_seqres"] = len(ent.full_sequence) if ent else None
    info["polymer_type"] = str(ent.polymer_type) if ent else None
    residues = []                                   # all polymer residues, first conformer
    for res in poly.first_conformer():
        ca = res.find_atom("CA", "*", gemmi.Element("C"))  # first altloc
        residues.append({
            "num": res.seqid.num, "icode": res.seqid.icode.strip(),
            "name": res.name, "label_seq": res.label_seq,
            "xyz": None if ca is None else (ca.pos.x, ca.pos.y, ca.pos.z)})
    return info, residues


def domain_nodes(segs_by_dom, residues, chain_char):
    """Map CATH segments (author numbering) to node ids. Returns list of node-id
    lists, or raises ValueError if a segment end is not found in the model."""
    key = {}
    for i, r in enumerate(residues):
        key.setdefault((r["num"], r["icode"]), i)
    node_of = {}
    for i, r in enumerate(residues):
        if r["xyz"] is not None:
            node_of[i] = len(node_of)
    groups = []
    for segs in segs_by_dom:
        nodes = []
        for c1, s, si, c2, e, ei in segs:
            if c1 != chain_char or c2 != chain_char:
                raise ValueError("segment on other chain")
            if (s, si) not in key or (e, ei) not in key:
                raise ValueError(f"segment end not in model: {s}{si}-{e}{ei}")
            a, b = key[(s, si)], key[(e, ei)]
            if b < a:
                raise ValueError("reversed segment")
            nodes.extend(node_of[i] for i in range(a, b + 1) if i in node_of)
        groups.append(sorted(set(nodes)))
    return groups


def contact_edges(xyz, cutoff):
    d2 = ((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1)
    iu, ju = np.where(np.triu(d2 <= cutoff * cutoff, k=1))
    return list(zip(iu.tolist(), ju.tolist()))


def stats(n, edges):
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    return {"m": g.number_of_edges(), "mean_degree": 2 * g.number_of_edges() / n,
            "mean_clustering": nx.average_clustering(g),
            "connected": nx.is_connected(g)}


# ----------------------------------------------------------------------- main
def main(raw, out):
    bounds = parse_boundaries(os.path.join(raw, "cath-domain-boundaries-v4_4_0.txt"))
    dlist = parse_domain_list(os.path.join(raw, "cath-domain-list-v4_4_0.txt"))
    cifdir = os.path.join(raw, "cif")
    os.makedirs(cifdir, exist_ok=True)
    os.makedirs(out, exist_ok=True)

    # Step 1: candidate chains from the CATH files alone.
    combos = defaultdict(list)
    n_multi = 0
    for chain, (doms, nf) in bounds.items():
        if len(doms) < 2:
            continue
        n_multi += 1
        if chain[4] == "0":            # no chain field in the PDB file
            continue
        ids = [f"{chain}{k + 1:02d}" for k in range(len(doms))]
        if not all(i in dlist for i in ids):
            continue
        res = dlist[ids[0]][2]
        if res > MAX_RES:              # also drops NMR (999) and obsolete (1000)
            continue
        tot = sum(dlist[i][1] for i in ids)
        if tot < MIN_LEN - 10 or tot > MAX_LEN:
            continue
        combo = tuple(sorted(dlist[i][0] for i in ids))
        combos[combo].append(chain)
    combo_keys = sorted(combos)
    rng = random.Random(SEED)
    rng.shuffle(combo_keys)
    print(f"multi-domain chains: {n_multi}; candidate combos: {len(combo_keys)}; "
          f"candidate chains: {sum(len(v) for v in combos.values())}", flush=True)

    # Step 2: walk combos in seeded random order; try up to 3 chains per combo.
    log = []
    chosen = []
    used_pdb = set()
    n_dl = 0
    for combo in combo_keys:
        if len(chosen) >= N_TARGET or n_dl >= MAX_DOWNLOADS:
            break
        chains = sorted(combos[combo])
        rng.shuffle(chains)
        tries = 0
        for chain in chains:
            if tries >= MAX_TRIES_PER_COMBO:
                break
            pdb, cc = chain[:4], chain[4]
            if pdb in used_pdb:
                log.append((chain, "|".join(combo), "skip: PDB entry already used"))
                continue
            tries += 1
            path, new = fetch(pdb, cifdir)
            n_dl += bool(new is True)
            if path is None:
                log.append((chain, "|".join(combo), f"reject: download failed ({new})"))
                continue
            try:
                info, residues = load_chain(path, cc)
            except Exception as exc:  # noqa: BLE001
                log.append((chain, "|".join(combo), f"reject: parse error {exc}"))
                continue
            reason = None
            if info["method"] != "X-RAY DIFFRACTION" or info["n_methods"] != 1:
                reason = f"method {info['method']}"
            elif not (info["resolution"] and info["resolution"] <= MAX_RES):
                reason = f"resolution {info['resolution']}"
            elif residues is None:
                reason = "chain not found"
            elif info.get("polymer_type") != "PolymerType.PeptideL":
                reason = f"polymer type {info.get('polymer_type')}"
            if reason is None:
                ca = [r for r in residues if r["xyz"] is not None]
                n = len(ca)
                labs = sorted(r["label_seq"] for r in ca if r["label_seq"] is not None)
                internal_missing = (labs[-1] - labs[0] + 1 - len(set(labs))) if labs else None
                total_missing = (info["n_seqres"] - len(set(labs))) if info["n_seqres"] else None
                if not (MIN_LEN <= n <= MAX_LEN):
                    reason = f"n_CA {n}"
                elif internal_missing is None or len(labs) != n:
                    reason = "label_seq missing"
                elif internal_missing > MAX_INTERNAL_MISSING:
                    reason = f"internal missing {internal_missing}"
                elif total_missing is None or total_missing > MAX_TOTAL_MISSING_FRAC * info["n_seqres"]:
                    reason = f"total missing {total_missing} of {info['n_seqres']}"
            if reason is None:
                try:
                    groups = domain_nodes(bounds[chain][0], residues, cc)
                except ValueError as exc:
                    reason = f"domain mapping: {exc}"
            if reason is None:
                kept = [g for g in groups if len(g) >= MIN_DOMAIN]
                allnodes = [v for g in kept for v in g]
                if len(allnodes) != len(set(allnodes)):
                    reason = "overlapping domains"
                elif len(kept) < 2:
                    reason = "fewer than 2 domains with >=3 residues"
                elif len(allnodes) < MIN_FRAC_IN_DOMAIN * n:
                    reason = f"domain coverage {len(allnodes) / n:.3f}"
                else:
                    xyz = np.array([r["xyz"] for r in ca], dtype=float)
                    if not all(stats(n, contact_edges(xyz, cut))["connected"] for cut in CUTOFFS):
                        reason = "contact network not connected"
            if reason is not None:
                log.append((chain, "|".join(combo), "reject: " + reason))
                continue
            log.append((chain, "|".join(combo), "accept"))
            used_pdb.add(pdb)
            chosen.append({
                "chain": chain, "combo": combo, "path": path, "info": info,
                "residues": [r for r in residues if r["xyz"] is not None],
                "groups": kept, "n_dom_cath": len(groups),
                "dom_ids": [f"{chain}{k + 1:02d}" for k in range(len(groups))],
                "dom_len_cath": [dlist[f"{chain}{k + 1:02d}"][1] for k in range(len(groups))],
                "internal_missing": internal_missing, "total_missing": total_missing,
                "n_frag": bounds[chain][1]})
            print(f"[{len(chosen):2d}] {chain} n={n} doms={[len(g) for g in kept]} "
                  f"res={info['resolution']} miss={internal_missing}/{total_missing}", flush=True)
            break

    # Step 3: write networks.
    rows = []
    for c in chosen:
        res = c["residues"]
        n = len(res)
        xyz = np.array([r["xyz"] for r in res], dtype=float)
        row = {"chain": c["chain"], "pdb": c["chain"][:4], "n": n,
               "resolution": c["info"]["resolution"],
               "n_domains": len(c["groups"]), "domain_sizes": [len(g) for g in c["groups"]],
               "cath_domain_lengths": c["dom_len_cath"],
               "superfamilies": list(c["combo"]),
               "frac_in_domain": sum(len(g) for g in c["groups"]) / n,
               "internal_missing": c["internal_missing"], "total_missing": c["total_missing"],
               "cath_fragments": c["n_frag"]}
        for cut in CUTOFFS:
            folder = c["chain"] if cut == 8.0 else f"{c['chain']}_7A"
            d = os.path.join(out, folder)
            os.makedirs(d, exist_ok=True)
            edges = contact_edges(xyz, cut)
            with open(os.path.join(d, "edges.tsv"), "w", newline="\n") as fh:
                fh.writelines(f"{u}\t{v}\n" for u, v in edges)
            with open(os.path.join(d, "groups.txt"), "w", newline="\n") as fh:
                fh.writelines(" ".join(map(str, g)) + "\n" for g in c["groups"])
            with open(os.path.join(d, "nodes.tsv"), "w", newline="\n") as fh:
                fh.write("id\tresidue_number\tresidue_name\n")
                fh.writelines(f"{i}\t{r['num']}{r['icode']}\t{r['name']}\n"
                              for i, r in enumerate(res))
            s = stats(n, edges)
            tag = "8A" if cut == 8.0 else "7A"
            for k, v in s.items():
                row[f"{k}_{tag}"] = v
        rows.append(row)

    with open(os.path.join(out, "selection_log.tsv"), "w", newline="\n") as fh:
        fh.write("chain\tsuperfamily_combo\toutcome\n")
        fh.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in log)
    summ = {"seed": SEED, "n_chains": len(rows), "downloads": n_dl,
            "candidate_combos": len(combo_keys),
            "candidate_chains": sum(len(v) for v in combos.values()),
            "chains": rows}
    for key in ("n", "m_8A", "mean_degree_8A", "mean_clustering_8A",
                "m_7A", "mean_degree_7A", "mean_clustering_7A", "n_domains", "frac_in_domain"):
        summ["mean_" + key] = float(np.mean([r[key] for r in rows])) if rows else None
    with open(os.path.join(out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "chains"}, indent=1))


def table(out):
    """Print the per-chain Markdown table for README.md from summary.json."""
    with open(os.path.join(out, "summary.json"), encoding="utf-8") as fh:
        s = json.load(fh)
    print("| chain | res. (A) | n | m (8 A) | domains | domain sizes | frac. in domain "
          "| mean deg. (8 A) | mean clust. (8 A) | m (7 A) | mean deg. (7 A) "
          "| mean clust. (7 A) | missing int./total | superfamilies |")
    print("|" + "---|" * 14)
    for r in s["chains"]:
        print(f"| {r['chain']} | {r['resolution']:.2f} | {r['n']} | {r['m_8A']} | {r['n_domains']} "
              f"| {' '.join(map(str, r['domain_sizes']))} | {r['frac_in_domain']:.3f} "
              f"| {r['mean_degree_8A']:.2f} | {r['mean_clustering_8A']:.3f} | {r['m_7A']} "
              f"| {r['mean_degree_7A']:.2f} | {r['mean_clustering_7A']:.3f} "
              f"| {r['internal_missing']}/{r['total_missing']} | {' '.join(r['superfamilies'])} |")
    print(f"| **mean over {s['n_chains']} chains** | |{s['mean_n']:.1f} | {s['mean_m_8A']:.1f} "
          f"| {s['mean_n_domains']:.2f} | | {s['mean_frac_in_domain']:.3f} "
          f"| {s['mean_mean_degree_8A']:.2f} | {s['mean_mean_clustering_8A']:.3f} "
          f"| {s['mean_m_7A']:.1f} | {s['mean_mean_degree_7A']:.2f} "
          f"| {s['mean_mean_clustering_7A']:.3f} | | |")


if __name__ == "__main__":
    if sys.argv[1] == "table":
        table(sys.argv[2])
    else:
        main(sys.argv[1], sys.argv[2])
