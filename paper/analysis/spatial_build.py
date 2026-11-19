"""Build spot graphs from 10x Visium spatial transcriptomics of human dorsolateral
prefrontal cortex (DLPFC), with the manual cortical-layer annotations as groups.

Source data (downloaded beforehand, plain formats only, never executed):
  * Maynard et al. 2021 (Nat. Neurosci. 24:425), 12 Visium slices 151507-151676.
    Space Ranger filtered count matrices (10x HDF5) from the LIBD public S3 bucket
    listed in the official LieberInstitute/HumanPilot repository
    (AWS_File_locations.tsv):
      https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/<S>_filtered_feature_bc_matrix.h5
  * Spot positions (Space Ranger tissue_positions_list.txt) and the barcode-level
    layer annotation (barcode_level_layer_map.tsv) from the same repository,
    pinned to commit 044446d6bd8fc154aa74f7be62ec67effb1ec376:
      https://raw.githubusercontent.com/LieberInstitute/HumanPilot/<sha>/10X/<S>/tissue_positions_list.txt
      https://raw.githubusercontent.com/LieberInstitute/HumanPilot/<sha>/10X/barcode_level_layer_map.tsv

RAW layout (each download in its own folder):
  RAW/h5_<S>/<S>_filtered_feature_bc_matrix.h5
  RAW/pos_<S>/tissue_positions_list.txt
  RAW/layer_map/barcode_level_layer_map.tsv

Usage (run from paper/, isolated interpreter; needs numpy, scipy, h5py, networkx,
scikit-learn for one diagnostic):
  python -I analysis/spatial_build.py inspect --raw D:/tmp/spatial_raw
  python -I analysis/spatial_build.py build   --raw D:/tmp/spatial_raw --out analysis/data_spatial
  (options: --samples ... --ks 8 16 --groups layers|parts|both --parts-ks 8)
  Default (--groups both --parts-ks 8): layer groups for every slice and k, written
  as dlpfc_<S>_k<K>[_null]; layer-part groups only for k = 8 and only on slices
  where the parts differ from the layers, written as dlpfc_<S>_parts_k8[_null].

Rule (fixed before any detection was run; see data_spatial/README.md):
  Nodes   spots in the filtered matrix that carry a layer label (L1-L6, WM),
          ordered by (array_row, array_col), ids 0..n-1.
  Feats   genes detected (count > 0) in >= 3 nodes; counts scaled to 10^4 per spot,
          log1p; 2000 highly variable genes (Seurat v1 dispersion rule, 20 mean bins,
          as scanpy flavor="seurat"); genes centred and scaled to unit variance,
          clipped at 10; PCA to 30 components (exact, eigendecomposition of Y'Y).
  Cands   for node i, the nodes at hexagonal grid distance 1, 2 or 3 (at most 36).
  Edges   i keeps its k candidates with the smallest Euclidean distance in PCA space
          (ties: smaller grid distance, then smaller id); undirected union.
  Null    i keeps its k candidates with the smallest Euclidean distance on the grid
          (ties: smaller id); undirected union. Same nodes, same groups.
  Groups  layers: one group per layer.  parts: each layer split into the connected
          components of the candidate graph (hex distance <= 3) restricted to that
          layer's nodes; parts with < 3 nodes are dropped.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # multithreaded BLAS was ~30x slower here
import numpy as np  # noqa: E402

SAMPLES = ["151507", "151508", "151509", "151510", "151669", "151670",
           "151671", "151672", "151673", "151674", "151675", "151676"]
LAYERS = ["L1", "L2", "L3", "L4", "L5", "L6", "WM"]
MIN_GENE_SPOTS = 3
TARGET_SUM = 1e4
N_HVG = 2000
N_PCS = 30
MAX_HEX = 3
MIN_GROUP = 3


# ----------------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------------
def read_layers(raw):
    out = {}
    with open(os.path.join(raw, "layer_map", "barcode_level_layer_map.tsv"),
              encoding="ascii") as fh:
        for line in fh:
            b, s, lay = line.rstrip("\n").split("\t")
            assert lay in LAYERS, lay
            assert (s, b) not in out
            out[(s, b)] = lay
    return out


def read_positions(raw, s):
    pos = {}
    with open(os.path.join(raw, f"pos_{s}", "tissue_positions_list.txt"),
              encoding="ascii") as fh:
        for row in csv.reader(fh):
            if row[0] == "barcode":
                continue
            b, it, r, c, pr, pc = row
            pos[b] = (int(it), int(r), int(c), int(pr), int(pc))
    return pos


def read_h5(raw, s):
    """Spots x genes CSR counts from a 10x v3 HDF5 file (stored genes x spots CSC)."""
    import h5py
    import scipy.sparse as sp
    with h5py.File(os.path.join(raw, f"h5_{s}", f"{s}_filtered_feature_bc_matrix.h5"),
                   "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][:])
        X = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=(shape[1], shape[0]))
        barcodes = [x.decode() for x in g["barcodes"][:]]
        genes = [x.decode() for x in g["features/id"][:]]
        ftype = [x.decode() for x in g["features/feature_type"][:]]
    assert set(ftype) == {"Gene Expression"}
    return X, barcodes, genes


def load_sample(raw, s, layers):
    X, barcodes, genes = read_h5(raw, s)
    pos = read_positions(raw, s)
    lab = {b: lay for (ss, b), lay in layers.items() if ss == s}
    bi = {b: i for i, b in enumerate(barcodes)}
    assert set(lab) <= set(bi), "labelled barcode missing from the matrix"
    assert all(pos[b][0] == 1 for b in barcodes), "matrix spot not in tissue"
    keep = sorted(lab, key=lambda b: (pos[b][1], pos[b][2]))
    rows = np.array([bi[b] for b in keep])
    return dict(
        sample=s, X=X[rows].astype(np.float64), barcodes=keep,
        layer=np.array([lab[b] for b in keep]),
        row=np.array([pos[b][1] for b in keep]),
        col=np.array([pos[b][2] for b in keep]),
        prow=np.array([pos[b][3] for b in keep]),
        pcol=np.array([pos[b][4] for b in keep]),
        n_matrix_spots=len(barcodes), n_in_tissue=sum(v[0] for v in pos.values()),
        n_genes_total=len(genes))


# ----------------------------------------------------------------------------
# expression features
# ----------------------------------------------------------------------------
def pca_features(X):
    import scipy.sparse as sp
    ngenes_det = np.asarray((X > 0).sum(0)).ravel()
    X = X[:, ngenes_det >= MIN_GENE_SPOTS]
    lib = np.asarray(X.sum(1)).ravel()
    assert (lib > 0).all()
    Xn = sp.diags(TARGET_SUM / lib) @ X
    Xl = Xn.log1p().tocsc()
    # HVG, Seurat v1 dispersion (scanpy flavor="seurat"): on expm1(log data)
    E = Xl.expm1()
    mean = np.asarray(E.mean(0)).ravel()
    var = np.asarray(E.multiply(E).mean(0)).ravel() - mean ** 2
    var *= E.shape[0] / (E.shape[0] - 1)
    mean[mean == 0] = 1e-12
    disp = var / mean
    disp[disp == 0] = np.nan
    disp = np.log(disp)
    lmean = np.log1p(mean)
    edges = np.linspace(lmean.min(), lmean.max(), 21)
    b = np.digitize(lmean, edges[1:-1], right=True)      # 20 equal-width bins, (a, b]
    nd = np.full_like(disp, np.nan)
    for k in range(20):
        m = b == k
        d = disp[m]
        ok = np.isfinite(d)
        if ok.sum() == 0:
            continue
        if ok.sum() == 1:                 # single-gene bin: scanpy sets mean 0, std = disp
            nd[m] = np.where(ok, 1.0, np.nan)
            continue
        nd[m] = (d - d[ok].mean()) / d[ok].std(ddof=1)
    order = np.argsort(np.where(np.isnan(nd), -np.inf, -nd), kind="stable")
    hvg = order[:N_HVG]
    Y = Xl[:, hvg].toarray()
    Y -= Y.mean(0)
    sd = Y.std(0, ddof=1)
    sd[sd == 0] = 1.0
    Y /= sd
    np.clip(Y, -np.inf, 10.0, out=Y)
    Y -= Y.mean(0)
    # exact PCA via the eigendecomposition of the 2000 x 2000 Gram matrix
    w, V = np.linalg.eigh(Y.T @ Y)
    w, V = w[::-1], V[:, ::-1]
    V = V[:, :N_PCS] * np.sign(V[np.abs(V[:, :N_PCS]).argmax(0), np.arange(N_PCS)])
    Z = Y @ V
    expl = float(w[:N_PCS].sum() / w.sum())
    return Z, dict(n_genes_filtered=int(X.shape[1]), n_hvg=len(hvg),
                   pca_var_explained=expl, median_umi=float(np.median(lib)))


# ----------------------------------------------------------------------------
# grid
# ----------------------------------------------------------------------------
def hex_dist(r1, c1, r2, c2):
    """Steps on the Visium hex grid (neighbours at (r, c+-2) and (r+-1, c+-1))."""
    dr = np.abs(r1 - r2)
    dc = np.abs(c1 - c2)
    return dr + np.maximum(0, (dc - dr) // 2)


def grid_euclid(r1, c1, r2, c2):
    """Centre-to-centre distance in units of the spot pitch."""
    return np.sqrt((0.5 * (c1 - c2)) ** 2 + (np.sqrt(3) / 2 * (r1 - r2)) ** 2)


def candidates(D):
    r, c = D["row"], D["col"]
    assert len(set((r + c) % 2)) == 1, "array_row + array_col parity not constant"
    n = len(r)
    index = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(r, c))}
    assert len(index) == n
    offs = []
    for dr in range(-MAX_HEX, MAX_HEX + 1):
        for dc in range(-2 * MAX_HEX, 2 * MAX_HEX + 1):
            if (dr + dc) % 2 or (dr == 0 and dc == 0):
                continue
            if hex_dist(0, 0, dr, dc) <= MAX_HEX:
                offs.append((dr, dc))
    assert len(offs) == 36
    cand = []
    for i in range(n):
        lst = [index[(int(r[i]) + dr, int(c[i]) + dc)] for dr, dc in offs
               if (int(r[i]) + dr, int(c[i]) + dc) in index]
        cand.append(np.array(sorted(lst), dtype=np.int64))
    return cand


def topk_edges(cand, keys_fn, k):
    E = set()
    for i, cj in enumerate(cand):
        if len(cj) == 0:
            continue
        keys = keys_fn(i, cj)          # tuple of arrays, primary key first
        order = np.lexsort(tuple(reversed(keys)) )
        for j in cj[order[:k]]:
            j = int(j)
            E.add((i, j) if i < j else (j, i))
    return sorted(E)


def layer_parts(D, cand):
    """Connected parts of each layer on the candidate graph restricted to that layer."""
    lay = D["layer"]
    n = len(lay)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in cand[i]:
            if lay[j] == lay[i]:
                a, b = find(i), find(int(j))
                if a != b:
                    parent[max(a, b)] = min(a, b)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), []).append(i)
    parts = sorted(comp.values(), key=lambda g: (LAYERS.index(lay[g[0]]), -len(g), g[0]))
    return parts


# ----------------------------------------------------------------------------
# stats
# ----------------------------------------------------------------------------
def frac_within(n, edges, groups):
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    return float(np.mean([bool(memb[a] & memb[b]) for a, b in edges]))


def graph_stats(n, edges, groups):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    sizes = [len(g) for g in groups]
    deg = np.array([d for _, d in G.degree()])
    return dict(n=n, m=len(edges), mean_degree=2 * len(edges) / n,
                min_degree=int(deg.min()), max_degree=int(deg.max()),
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                n_isolated=int((deg == 0).sum()),
                n_groups=len(groups), group_size_min=min(sizes),
                group_size_median=float(np.median(sizes)), group_size_max=max(sizes),
                frac_nodes_in_group=len(set().union(*map(set, groups))) / n,
                frac_edges_within_group=float(np.mean(
                    [bool(memb[a] & memb[b]) for a, b in edges])))


def write_folder(od, edges, groups, D):
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "edges.tsv"), "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)
    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
        for g in groups:
            fh.write(" ".join(map(str, sorted(g))) + "\n")
    with open(os.path.join(od, "nodes.tsv"), "w", newline="\n") as fh:
        fh.write("id\tbarcode\tarray_row\tarray_col\tpxl_row_in_fullres\t"
                 "pxl_col_in_fullres\tlayer\n")
        for i, b in enumerate(D["barcodes"]):
            fh.write(f"{i}\t{b}\t{D['row'][i]}\t{D['col'][i]}\t{D['prow'][i]}\t"
                     f"{D['pcol'][i]}\t{D['layer'][i]}\n")


# ----------------------------------------------------------------------------
# commands
# ----------------------------------------------------------------------------
def cmd_inspect(args):
    layers = read_layers(args.raw)
    for s in args.samples:
        D = load_sample(args.raw, s, layers)
        cand = candidates(D)
        parts = layer_parts(D, cand)
        nc = np.array([len(c) for c in cand])
        by = {}
        for p in parts:
            by.setdefault(D["layer"][p[0]], []).append(len(p))
        print(s, "n", len(D["barcodes"]), "matrix", D["n_matrix_spots"],
              "cands min/median", nc.min(), np.median(nc),
              {k: v[:6] + (["..."] if len(v) > 6 else []) for k, v in by.items()},
              "parts>=3:", sum(len(p) >= MIN_GROUP for p in parts),
              "nodes in parts>=3:", sum(len(p) for p in parts if len(p) >= MIN_GROUP))


def cmd_build(args):
    layers = read_layers(args.raw)
    os.makedirs(args.out, exist_ok=True)
    summary = []
    for s in args.samples:
        D = load_sample(args.raw, s, layers)
        n = len(D["barcodes"])
        Z, finfo = pca_features(D["X"])
        cand = candidates(D)
        r, c = D["row"], D["col"]
        group_sets = {}
        if args.groups in ("layers", "both"):
            group_sets["layers"] = [list(np.nonzero(D["layer"] == L)[0])
                                    for L in LAYERS if (D["layer"] == L).sum() >= MIN_GROUP]
        parts_differ = None
        if args.groups in ("parts", "both"):
            parts = [p for p in layer_parts(D, cand) if len(p) >= MIN_GROUP]
            lay = [list(np.nonzero(D["layer"] == L)[0]) for L in LAYERS
                   if (D["layer"] == L).sum() >= MIN_GROUP]
            parts_differ = sorted(map(tuple, parts)) != sorted(map(tuple, lay))
            # a parts set is written only where it differs from the layers
            if parts_differ or args.groups == "parts":
                group_sets["parts"] = parts

        def expr_keys(i, cj):
            return (np.linalg.norm(Z[cj] - Z[i], axis=1),
                    grid_euclid(r[i], c[i], r[cj], c[cj]), cj)

        def grid_keys(i, cj):
            return (grid_euclid(r[i], c[i], r[cj], c[cj]), cj)

        rng = np.random.default_rng(int(s))
        rscore = [rng.random(len(x)) for x in cand]

        def rand_keys(i, cj):              # diagnostic only, never written as a graph
            return (rscore[i], cj)

        nc = np.array([len(x) for x in cand])
        cand_pairs = [(i, int(j)) for i in range(n) for j in cand[i] if i < j]
        for k in args.ks:
            E = topk_edges(cand, expr_keys, k)
            N = topk_edges(cand, grid_keys, k)
            Rn = topk_edges(cand, rand_keys, k)
            # diagnostic only: expression kNN with no spatial restriction
            from sklearn.neighbors import NearestNeighbors
            _, nb = NearestNeighbors(n_neighbors=k + 1).fit(Z).kneighbors(Z)
            G_all = sorted({(min(i, int(j)), max(i, int(j)))
                            for i in range(n) for j in nb[i] if int(j) != i})
            hd = hex_dist(r[[a for a, _ in E]], c[[a for a, _ in E]],
                          r[[b for _, b in E]], c[[b for _, b in E]])
            for gname, groups in group_sets.items():
                if gname == "parts" and args.groups == "both" and k not in args.parts_ks:
                    continue
                tag = f"dlpfc_{s}_k{k}" if gname == "layers" else f"dlpfc_{s}_parts_k{k}"
                write_folder(os.path.join(args.out, tag), E, groups, D)
                write_folder(os.path.join(args.out, tag + "_null"), N, groups, D)
                for name, edges, ctl in ((tag, E, False), (tag + "_null", N, True)):
                    st = graph_stats(n, edges, groups)
                    row = dict(graph=name, sample=s, k=k, groups=gname, control=ctl,
                               rule="grid" if ctl else "expression", **st)
                    if not ctl:
                        nst = graph_stats(n, N, groups)
                        row.update(
                            null_m=nst["m"], null_clustering=nst["clustering"],
                            null_frac_edges_within_group=nst["frac_edges_within_group"],
                            null_n_components=nst["n_components"],
                            frac_edges_hex1=float(np.mean(hd == 1)),
                            frac_edges_hex2=float(np.mean(hd == 2)),
                            frac_edges_hex3=float(np.mean(hd == 3)),
                            frac_edges_shared_with_null=len(set(E) & set(N)) / len(E),
                            random_rule_m=len(Rn),
                            random_rule_frac_edges_within_group=frac_within(n, Rn, groups),
                            global_expr_knn_m=len(G_all),
                            global_expr_knn_frac_edges_within_group=frac_within(
                                n, G_all, groups),
                            candidate_pairs=len(cand_pairs),
                            candidate_pairs_frac_within_group=frac_within(
                                n, cand_pairs, groups),
                            n_spots_in_tissue=D["n_in_tissue"],
                            n_spots_matrix=D["n_matrix_spots"],
                            n_unlabelled_dropped=D["n_matrix_spots"] - n,
                            n_layers=int(len(set(D["layer"]))),
                            layer_sizes={L: int((D["layer"] == L).sum())
                                         for L in LAYERS if (D["layer"] == L).any()},
                            candidates_min=int(nc.min()),
                            candidates_median=float(np.median(nc)),
                            n_nodes_lt_k_candidates=int((nc < k).sum()),
                            layer_parts_differ=parts_differ, **finfo)
                    summary.append(row)
                    print(json.dumps({x: row[x] for x in (
                        "graph", "n", "m", "mean_degree", "clustering", "n_components",
                        "n_groups", "frac_edges_within_group")}), flush=True)
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("inspect", "build"):
        a = sub.add_parser(name)
        a.add_argument("--raw", required=True)
        a.add_argument("--samples", nargs="+", default=SAMPLES)
        if name == "build":
            a.add_argument("--out", required=True)
            a.add_argument("--ks", nargs="+", type=int, default=[8, 16])
            a.add_argument("--groups", choices=["layers", "parts", "both"],
                           default="both")
            a.add_argument("--parts-ks", nargs="+", type=int, default=[8],
                           help="with --groups both: k values that also get a parts set")
    args = ap.parse_args()
    {"inspect": cmd_inspect, "build": cmd_build}[args.cmd](args)


if __name__ == "__main__":
    main()
