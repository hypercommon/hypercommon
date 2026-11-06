"""
Experiment definitions. Each family is a fixed graph model; p is the rewired
fraction. Fixed here, before results, so the design is predeclared.
"""

import graphs

ALL_METHODS = ["hc", "hc_2test", "hc_tri", "cpm", "leiden", "infomap", "lpa",
               "walktrap", "slpa", "demon", "angel", "scan", "scan_j", "linkcomm"]
P_GRID = [round(0.05 * i, 2) for i in range(17)]          # 0 .. 0.80


def _rings(fam, p, seed):
    return graphs.rings(fam["sizes"], fam["ks"], fam["overlap"], p, seed,
                        rewire=fam.get("rewire", "uniform"), internal=fam.get("internal", 0.0),
                        isolated=fam.get("isolated", False),
                        internal_mode=fam.get("internal_mode", "endpoint"),
                        split=fam.get("split", False), hubs=fam.get("hubs"),
                        geometry=fam.get("geometry", "lattice"))


def _null(fam, p, seed):
    return graphs.null_graph(fam["kind"], seed)


def _uniform(n_rings, size, k, overlap):
    return {"sizes": [size] * n_rings, "ks": [k] * n_rings, "overlap": overlap}


def _shape(parts, overlap):
    sizes, ks = [], []
    for size, k, count in parts:
        sizes += [size] * count
        ks += [k] * count
    return {"sizes": sizes, "ks": ks, "overlap": overlap}


def _lfr(fam, mu, seed):
    return graphs.lfr(fam, mu, seed)


LFR_BASE = dict(N=2000, k=20, maxk=50, t1=2, t2=1, minc=20, maxc=100, om=2)
MU_GRID = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

EXPERIMENTS = {
    # Uniform rings: the setting where modularity methods are strongest.
    "rings_uniform": {
        "build": _rings,
        "families": {
            f"k{k}_s{s}_o{int(o*100):02d}": _uniform(2000 // s, s, k, o)
            for (k, s) in [(8, 50), (16, 100)]
            for o in (0.0, 0.05, 0.10)
        },
        "p_grid": P_GRID,
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # Heterogeneous rings: sizes and degrees both vary (from shapes_* sweeps),
    # including one dense, even shape where Leiden previously won.
    "rings_hetero": {
        "build": _rings,
        "families": {
            "two_worlds": _shape([(200, 32, 5), (50, 8, 20)], 0.05),
            "three_big_many_small": _shape([(400, 32, 3), (50, 8, 16)], 0.05),
            "gradient10": _shape([(600, 16, 1), (400, 12, 1), (300, 12, 1), (200, 8, 1), (150, 8, 1),
                                  (120, 8, 1), (100, 8, 1), (60, 8, 1), (40, 8, 1), (30, 8, 1)], 0.05),
            "dense_even": _shape([(100, 16, 10), (50, 16, 20)], 0.05),     # same k, sizes differ
            "sizes_only_k8": _shape([(100, 8, 10), (50, 8, 20)], 0.05),    # same k, sizes differ
            "k_only_s100": _shape([(100, 32, 5), (100, 24, 5), (100, 16, 5), (100, 8, 5)], 0.05),  # same s, k differ
        },
        "p_grid": P_GRID,
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # Overlapping LFR (Lancichinetti & Fortunato 2009); the noise axis is mu.
    "lfr": {
        "build": _lfr,
        "families": {f"on{int(f*100):02d}": dict(LFR_BASE, on_frac=f) for f in (0.0, 0.1, 0.3)},
        "p_grid": MU_GRID,
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # One factor at a time around (mu, on=0.1): sparser, more memberships, larger communities.
    "lfr_slices": {
        "build": _lfr,
        "families": {
            "k10": dict(LFR_BASE, k=10, maxk=25, on_frac=0.1),
            "om4": dict(LFR_BASE, om=4, on_frac=0.1),
            "big": dict(LFR_BASE, minc=50, maxc=250, on_frac=0.1),
        },
        "p_grid": [0.2, 0.4, 0.6],
        "p_grid_v2": MU_GRID,               # main runs: the same mu range as the other LFR families
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # Commonality as a FAMILY: HyperCommon with each f (quantile threshold grid),
    # Omega only. Graph families reuse the definitions above.
    "commonality_family": {
        "build": lambda fam, p, seed: (_rings if fam["kind"] == "rings" else _lfr)(fam, p, seed),
        "families": {
            "rings_k8_o05": dict(_uniform(40, 50, 8, 0.05), kind="rings"),
            "rings_k16_o05": dict(_uniform(20, 100, 16, 0.05), kind="rings"),
            "two_worlds": dict(_shape([(200, 32, 5), (50, 8, 20)], 0.05), kind="rings"),
            "gradient10": dict(_shape([(600, 16, 1), (400, 12, 1), (300, 12, 1), (200, 8, 1), (150, 8, 1),
                                       (120, 8, 1), (100, 8, 1), (60, 8, 1), (40, 8, 1), (30, 8, 1)], 0.05), kind="rings"),
            "lfr_on10": dict(LFR_BASE, on_frac=0.1, kind="lfr"),
            "lfr_on30": dict(LFR_BASE, on_frac=0.3, kind="lfr"),
            "lfr_k10": dict(LFR_BASE, k=10, maxk=25, on_frac=0.1, kind="lfr"),
        },
        "p_grid": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        "p_grid_by_family": {"lfr_on10": MU_GRID, "lfr_on30": MU_GRID, "lfr_k10": MU_GRID},
        "methods": ["hcf_" + f for f in ("jaccard", "count", "i2_over_u", "i_over_sqrtu",
                                          "cosine", "min_norm", "max_norm", "resource_alloc")],
        "onmi": False,
    },
    # Degree-biased (preferential) rewiring: the moved endpoint attaches with
    # probability proportional to deg+1, so noise concentrates on hubs.
    # Unequal rings at overlap 0 and 0.10 (rings_hetero uses 0.05), so overlap is
    # varied on the families where size and density differ.
    "rings_hetero_overlap": {
        "build": _rings,
        "families": {
            f"{name}_o{int(o * 100):02d}": _shape(parts, o)
            for name, parts in (("two_worlds", [(200, 32, 5), (50, 8, 20)]),
                                ("three_big_many_small", [(400, 32, 3), (50, 8, 16)]),
                                ("gradient10", [(600, 16, 1), (400, 12, 1), (300, 12, 1), (200, 8, 1), (150, 8, 1),
                                                (120, 8, 1), (100, 8, 1), (60, 8, 1), (40, 8, 1), (30, 8, 1)]),
                                ("sizes_only_k8", [(100, 8, 10), (50, 8, 20)]),
                                # the two remaining unequal families
                                ("dense_even", [(100, 16, 10), (50, 16, 20)]),
                                ("k_only_s100", [(100, 32, 5), (100, 24, 5), (100, 16, 5), (100, 8, 5)]))
            for o in (0.0, 0.10)
        },
        "p_grid": P_GRID,
        "methods": ALL_METHODS + ["slpa_x", "demon_x", "angel_x", "cpm_x", "scan_x"],
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa", "angel_x", "demon_x", "slpa_x"],
    },
    "rings_pref": {
        "build": _rings,
        "families": {
            f"{name}_{rw}": dict(fam, rewire=rw)
            for name, fam in (("k8_s50_o05", _uniform(40, 50, 8, 0.05)),
                              ("two_worlds", _shape([(200, 32, 5), (50, 8, 20)], 0.05)))
            for rw in ("preferential", "preferential2")
        },
        "p_grid": [round(0.05 * i, 2) for i in range(21)],          # 0 .. 1.0
        "p_grid_v2": P_GRID,                                          # main runs: 0 .. 0.8
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # Community GEOMETRY at fixed sizes, degrees and overlap: a fraction q of each
    # ring's edges is first rewired inside the ring, from the lattice (q = 0) through
    # small-world, clustered communities to random graphs (q = 1); then merges and
    # external noise p as usual. Separates "unequal sizes" from "lattice shape".
    "rings_geometry": {
        "build": _rings,
        "families": {
            f"{name}_q{int(q * 100):03d}": dict(fam, internal=q)
            for name, fam in (("k8_s50_o05", _uniform(40, 50, 8, 0.05)),
                              ("two_worlds", _shape([(200, 32, 5), (50, 8, 20)], 0.05)),
                              ("three_big_many_small", _shape([(400, 32, 3), (50, 8, 16)], 0.05)),
                              ("gradient10", _shape([(600, 16, 1), (400, 12, 1), (300, 12, 1), (200, 8, 1),
                                                     (150, 8, 1), (120, 8, 1), (100, 8, 1), (60, 8, 1),
                                                     (40, 8, 1), (30, 8, 1)], 0.05)))
            for q in (0.0, 0.1, 0.3, 0.6, 1.0)
        },
        "p_grid": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        "methods": ALL_METHODS,
        "slow_first": ["walktrap", "angel", "demon", "slpa"],
    },
    # Theorem 2's hypothesis made true: merges more than 2k apart on every ring
    # (possible only at low overlap), to test the predicted exact-recovery window.
    "rings_isolated": {
        "build": _rings,
        "families": {f"k{k}_s{s}_o01_iso": dict(_uniform(2000 // s, s, k, 0.01), isolated=True)
                     for (k, s) in [(8, 50), (16, 100)]},
        "p_grid": [0.0, 0.05, 0.1, 0.2],
        "methods": ["hc", "hc_tri", "hc_2test"],
    },
    # LFR with realistic clustering: standard LFR loses its triangles as mu grows
    # (avg clustering 0.37 at mu=0.1, 0.07 at mu=0.5); -C holds it fixed.
    "lfr_clustered": {
        "build": _lfr,
        "families": {f"on10_C{int(C*10):02d}": dict(LFR_BASE, on_frac=0.1, C=C) for C in (0.3, 0.5)},
        "p_grid": MU_GRID,
        "methods": ALL_METHODS,
        "slow_first": ["linkcomm", "walktrap", "angel", "demon", "slpa"],
    },
    # ---- control and null families ------------------------------------------------
    # Hubs INSIDE clustered communities: ring lattices whose nodes draw their own
    # degree from a power law (gamma = 2.5, k_r .. 8 k_r), so communities keep their
    # lattice neighbourhoods but contain hubs (graphs._hub_lattice).
    "rings_hubs": {
        "build": _rings,
        "families": {
            "hub_s50": dict(_uniform(40, 50, 4, 0.05), hubs=(2.5, 8)),
            "hub_s100": dict(_uniform(20, 100, 8, 0.05), hubs=(2.5, 8)),
            "hub_two_worlds": dict(_shape([(200, 16, 5), (50, 4, 20)], 0.05), hubs=(2.5, 8)),
        },
        "p_grid": [round(0.1 * i, 1) for i in range(9)],
        "methods": ALL_METHODS,
    },
    # Degree-preserving control of the geometry families (double-edge swaps inside rings).
    "rings_geometry_swap": {
        "build": _rings,
        "families": {
            f"{name}_q{int(q * 100):03d}": dict(fam, internal=q, internal_mode="swap")
            for name, fam in (("k8_s50_o05", _uniform(40, 50, 8, 0.05)),
                              ("three_big_many_small", _shape([(400, 32, 3), (50, 8, 16)], 0.05)))
            for q in (0.1, 0.3)
        },
        "p_grid": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        "methods": ALL_METHODS,
    },
    # Shared nodes that split their degree between the two rings (as LFR's do).
    "rings_split": {
        "build": _rings,
        "families": {
            f"{name}_o{int(o * 100):02d}_split": dict(_shape(parts, o), split=True)
            for name, parts in (("two_worlds", [(200, 32, 5), (50, 8, 20)]),
                                ("three_big_many_small", [(400, 32, 3), (50, 8, 16)]))
            for o in (0.05, 0.10)
        },
        "p_grid": [round(0.1 * i, 1) for i in range(9)],
        "methods": ALL_METHODS,
    },
    # A clustered benchmark that is not a ring: random geometric communities.
    "rgg": {
        "build": _rings,
        "families": {
            "rgg_uniform": dict(_uniform(40, 50, 8, 0.05), geometry="rgg"),
            "rgg_two_worlds": dict(_shape([(200, 32, 5), (50, 8, 20)], 0.05), geometry="rgg"),
        },
        "p_grid": [round(0.1 * i, 1) for i in range(9)],
        "methods": ALL_METHODS,
    },
    # Clustering without communities: what each method reports when there is nothing.
    "nulls": {
        "build": _null,
        "families": {"rgg_null": {"kind": "rgg"}, "hk_null": {"kind": "holme_kim"}},
        "p_grid": [0.0],
        "methods": ALL_METHODS,
    },
}

