# Protein residue contact networks with CATH domains as ground truth

60 protein chains. Each chain is a graph: nodes are residues, edges are
C-alpha contacts. The ground-truth groups are the chain's CATH structural
domains. Built by `../protein_build.py` (run with `python -I`).

## Sources

| what | URL | version |
|---|---|---|
| CATH domain boundaries (CathDomall format 2.0) | http://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-boundaries-v4_4_0.txt | CATH v4.4.0, file date 16.12.2024, sha256 `0505e9dc...ca34da` |
| CATH domain list (superfamily code, domain length, resolution) | http://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-v4_4_0.txt | CATH v4.4.0, file date 16.12.2024, sha256 `2f87544b...e088cc81` |
| File-format notes | same folder: `README-cath-domain-boundaries-file-format.txt`, `README-cath-list-file-format.txt` | v4.4.0 |
| Atomic coordinates (PDBx/mmCIF) | `https://files.rcsb.org/download/<PDBID>.cif.gz` | RCSB PDB archive as downloaded on 2026-11-12 |

v4.4.0 is the newest numbered release on the CATH server
(`latest-release/` holds the same files without the version suffix).

Licence / terms:
- CATH: Creative Commons Attribution 4.0 International (CC BY 4.0), as stated on
  https://www.cathdb.info/about. Cite CATH (Sillitoe et al.) when using these groups.
- PDB coordinates: CC0 1.0 public domain dedication (wwPDB policy,
  https://www.rcsb.org/pages/usage-policy). Cite the individual PDB entries
  if they are discussed.
- No account or login was used. Downloads were sequential with at least 1 s
  between requests (130 mmCIF files: 129 fetched by the script, 1 manual test).

## Files

For each chain, folder `<chain>/` (8 A edges) and `<chain>_7A/` (7 A edges):

- `edges.tsv`: one undirected edge per line, `u<TAB>v`, with `u < v`; node ids
  0..n-1; no self-loops, no duplicates.
- `groups.txt`: one line per CATH domain (in CATH domain-number order),
  space-separated node ids. Discontinuous domains include all their segments.
  Residues outside every domain (CATH "fragments", linkers, unchopped tails)
  are in no group.
- `nodes.tsv`: header `id residue_number residue_name`; residue number is the
  author number (`auth_seq_id`) followed by the insertion code, if any (e.g. `132E`).

`groups.txt` and `nodes.tsv` are the same in both folders of a chain.
Also: `summary.json` (all per-chain numbers below, plus means) and
`selection_log.tsv` (every chain that was looked at, and why it was kept or
rejected). The chain id is the CATH chain id: 4-character PDB id + author chain id.

## Nodes and edges

- Structure: first model only. Residues: the polymer (L-peptide entity) of the
  author chain, first conformer only (for microheterogeneity, the first
  residue listed). Ligands and waters are excluded.
- Node: every polymer residue that has a C-alpha atom (atom name `CA`, element
  C), in chain order. Modified residues inside the chain (e.g. MSE) are kept.
  For alternate locations, the first `CA` atom listed (normally altloc A) is used.
- Edge rule: residues i != j are joined when the distance between their
  C-alpha atoms is <= 8.0 A (main version) or <= 7.0 A (robustness version,
  `_7A` folders). No sequence-separation exclusion: consecutive residues
  (C-alpha distance ~3.8 A) and i,i+2 neighbours are always joined.
- Domains: each CATH segment `start-end` (author number + insertion code) is
  located by exact match in the model's residue list; all residues from start
  to end in chain order, that have a C-alpha atom, go into that domain. Domains
  with fewer than 3 such residues would be dropped (none were).
- CATH's own domain lengths exclude HETATM residues such as MSE, so our domain
  sizes can be a few residues larger than the CATH length (at most 8 here).

## Selection

Seed: `2026` (Python `random.Random(2026)`). Target: 60 chains.

1. From the CATH v4.4.0 boundaries file, take chains with >= 2 domains
   (179,099 chains), skip chains with no chain id (`0`), and keep those whose
   CATH resolution is <= 2.5 A (this also drops NMR = 999 and obsolete = 1000)
   and whose summed CATH domain length is 140-800. This gives 81,052 chains.
2. Group these chains by superfamily combination: the sorted list (with
   repeats) of the CATH superfamily codes of the chain's domains. 2,830
   combinations.
3. Sort the combinations, shuffle them with the seeded RNG, and walk through
   them. For each combination, sort its chains, shuffle them with the same
   RNG, and try at most 3 chains (chains from an already used PDB entry are
   skipped without counting as a try). Take the first chain that passes all
   checks below, then move to the next combination. Stop at 60 chains.
4. Checks on the downloaded mmCIF:
   - `_exptl.method` is exactly one method, `X-RAY DIFFRACTION`, and the
     mmCIF resolution is <= 2.5 A (so the filter was available and applied);
   - 150 <= n <= 800, where n = residues with a C-alpha atom;
   - at most 5 residues missing inside the chain (SEQRES positions between
     the first and last observed C-alpha that have no C-alpha);
   - at most 10% of the SEQRES length unobserved in total (termini included);
   - all CATH segment ends are found in the model, domains do not overlap,
     and >= 2 domains have >= 3 residues;
   - at least 80% of the n residues lie in a domain;
   - the contact network is connected at both 7 A and 8 A.
5. Result: one chain per superfamily combination and at most one chain per
   PDB entry. 147 chains (129 PDB entries) were tried in 94 combinations (see
   `selection_log.tsv`): 60 accepted; rejected for internal gaps 63, total
   missing 9, length 6, domain coverage 4, not X-ray (3 EM, 1 joint
   X-ray/neutron) 4, disconnected network 1 (2istA: its N-terminal domain has
   no C-alpha within 8 A of the rest of the chain).

Selected chains (in selection order):

```
1hyuA 1c97A 3cloA 2wqyN 3earA 1s20F 4yy8A 4zi8A 1aisB 1ybaC 1zqwA 3ib3B 
3rgcA 1sg4C 5t3oC 3fjyB 1nbwA 3b55A 1e8gA 4pqaA 5axiB 5nswC 5i3zD 2e8gA 
2jjlA 1ro2A 2x49A 1h3jB 7k76L 2a7kA 3ihmB 1rysB 3cw7D 1ulyA 6v06A 1gz0A 
1b7yB 3fulA 6u45A 4ihfA 3sluB 2bnoA 4c73B 1qupA 2d1cA 2zo9B 1nxzB 1cyyA 
1v6vB 3smzA 1xdpA 3mogA 5le8B 1x0cB 3m7nB 2nv9F 5z24A 1it8A 4yb5D 2fjrB 
```

## Rebuild

```
python -I paper/analysis/protein_build.py RAW_DIR paper/analysis/data_protein
python -I paper/analysis/protein_build.py table paper/analysis/data_protein   # per-chain table
```

`RAW_DIR` must hold the two CATH files above; mmCIF files are fetched into
`RAW_DIR/cif/` (cached). Built with Python 3.12.10, gemmi 0.7.5, numpy 2.5.4,
networkx 3.7. The RCSB archive is updated weekly, so a later rebuild could in
principle see remediated coordinates; the selection itself depends only on the
CATH files, the seed and the checks.

## Summary (mean over the 60 chains)

| quantity | 8 A (main) | 7 A |
|---|---|---|
| n (nodes) | 369.4 | 369.4 |
| m (edges) | 1852.9 | 1497.7 |
| mean degree | 9.92 | 8.05 |
| mean clustering (networkx `average_clustering`) | 0.545 | 0.549 |
| number of domains | 2.78 | 2.78 |
| fraction of residues in a domain | 0.983 | 0.983 |

n ranges 179-775 (median 333); domains per chain: 2 (30 chains), 3 (17),
4 (10), 5 (2), 6 (1); domain sizes 25-384 residues; resolution 1.30-2.50 A.
All 120 graphs are connected. 8 chains have all domains in the same
superfamily (internal repeats or duplications, e.g. two Ig domains).

## Per-chain table

Mean degree = 2m/n. Mean clustering = mean local clustering coefficient
(nodes of degree < 2 count as 0). "frac. in domain" = fraction of the n
residues in a kept domain. "missing int./total" = SEQRES residues without
C-alpha inside the chain / in total including termini.

| chain | res. (A) | n | m (8 A) | domains | domain sizes | frac. in domain | mean deg. (8 A) | mean clust. (8 A) | m (7 A) | mean deg. (7 A) | mean clust. (7 A) | missing int./total | superfamilies |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1hyuA | 2.00 | 521 | 2639 | 4 | 182 95 123 94 | 0.948 | 10.13 | 0.534 | 2141 | 8.22 | 0.536 | 0/0 | 3.40.30.10 3.40.30.10 3.50.50.60 3.50.50.60 |
| 1c97A | 1.98 | 753 | 4098 | 4 | 201 113 175 221 | 0.943 | 10.88 | 0.514 | 3241 | 8.61 | 0.517 | 0/0 | 3.20.19.10 3.30.499.10 3.30.499.10 3.40.1060.10 |
| 3cloA | 2.04 | 258 | 1221 | 2 | 194 64 | 1.000 | 9.47 | 0.544 | 1016 | 7.88 | 0.562 | 0/0 | 1.10.10.10 3.30.450.20 |
| 2wqyN | 2.10 | 612 | 3300 | 4 | 326 104 117 65 | 1.000 | 10.78 | 0.522 | 2610 | 8.53 | 0.522 | 0/9 | 1.20.58.100 3.50.50.60 3.90.700.10 4.10.80.40 |
| 3earA | 2.30 | 206 | 956 | 2 | 153 53 | 1.000 | 9.28 | 0.575 | 813 | 7.89 | 0.590 | 3/6 | 3.40.50.300 6.10.250.710 |
| 1s20F | 2.20 | 333 | 1741 | 3 | 80 149 40 | 0.808 | 10.46 | 0.555 | 1417 | 8.51 | 0.553 | 0/7 | 1.10.1530.10 3.30.1370.60 3.30.60.50 |
| 4yy8A | 1.81 | 373 | 1906 | 2 | 75 298 | 1.000 | 10.22 | 0.554 | 1544 | 8.28 | 0.549 | 4/16 | 2.120.10.80 3.30.710.10 |
| 4zi8A | 1.70 | 313 | 1507 | 3 | 97 108 108 | 1.000 | 9.63 | 0.552 | 1224 | 7.82 | 0.547 | 0/21 | 2.60.40.60 2.60.40.60 2.60.40.60 |
| 1aisB | 2.10 | 193 | 908 | 2 | 99 94 | 1.000 | 9.41 | 0.563 | 768 | 7.96 | 0.581 | 0/7 | 1.10.472.10 1.10.472.10 |
| 1ybaC | 2.24 | 406 | 2008 | 3 | 133 189 82 | 0.995 | 9.89 | 0.537 | 1619 | 7.98 | 0.538 | 0/4 | 3.30.70.260 3.40.50.720 3.40.50.720 |
| 1zqwA | 2.30 | 242 | 1129 | 3 | 58 114 69 | 0.996 | 9.33 | 0.568 | 935 | 7.73 | 0.570 | 3/6 | 1.10.150.20 3.30.210.10 3.30.460.10 |
| 3ib3B | 2.05 | 552 | 2935 | 3 | 284 52 180 | 0.935 | 10.63 | 0.519 | 2316 | 8.39 | 0.511 | 4/8 | 1.10.3020.20 2.60.120.260 3.40.50.1820 |
| 3rgcA | 2.30 | 248 | 1111 | 2 | 170 78 | 1.000 | 8.96 | 0.588 | 960 | 7.74 | 0.592 | 3/4 | 1.10.4030.10 3.10.50.40 |
| 1sg4C | 1.30 | 250 | 1256 | 2 | 202 48 | 1.000 | 10.05 | 0.551 | 1051 | 8.41 | 0.564 | 0/10 | 3.90.226.10 6.10.250.170 |
| 5t3oC | 2.20 | 307 | 1507 | 2 | 165 139 | 0.990 | 9.82 | 0.542 | 1200 | 7.82 | 0.545 | 0/0 | 3.40.50.2020 3.40.50.2020 |
| 3fjyB | 2.15 | 352 | 1786 | 2 | 180 172 | 1.000 | 10.15 | 0.530 | 1457 | 8.28 | 0.543 | 0/12 | 3.40.50.1240 3.90.79.10 |
| 1nbwA | 2.40 | 606 | 3172 | 5 | 132 143 154 132 45 | 1.000 | 10.47 | 0.538 | 2531 | 8.35 | 0.535 | 0/1 | 2.40.50.140 3.30.420.40 3.30.420.40 3.50.30.70 3.90.470.30 |
| 3b55A | 2.30 | 408 | 2044 | 3 | 209 72 127 | 1.000 | 10.02 | 0.533 | 1660 | 8.14 | 0.541 | 0/43 | 1.20.1440.30 3.30.1870.10 3.40.1660.10 |
| 1e8gA | 2.10 | 550 | 2807 | 4 | 118 111 231 53 | 0.933 | 10.21 | 0.527 | 2245 | 8.16 | 0.534 | 5/10 | 1.10.45.10 3.30.43.10 3.30.465.10 3.40.462.10 |
| 4pqaA | 1.78 | 375 | 1938 | 2 | 261 114 | 1.000 | 10.34 | 0.534 | 1554 | 8.29 | 0.533 | 0/6 | 3.30.70.360 3.40.630.10 |
| 5axiB | 2.50 | 302 | 1417 | 3 | 128 89 85 | 1.000 | 9.38 | 0.557 | 1181 | 7.82 | 0.581 | 0/5 | 1.10.238.10 1.20.930.20 3.30.505.10 |
| 5nswC | 2.50 | 451 | 2169 | 2 | 370 81 | 1.000 | 9.62 | 0.549 | 1855 | 8.23 | 0.559 | 0/2 | 1.20.1600.10 2.20.200.10 |
| 5i3zD | 2.05 | 325 | 1676 | 2 | 216 109 | 1.000 | 10.31 | 0.517 | 1312 | 8.07 | 0.526 | 0/3 | 3.40.50.1170 3.40.50.40 |
| 2e8gA | 1.70 | 240 | 1181 | 2 | 124 116 | 1.000 | 9.84 | 0.543 | 963 | 8.03 | 0.546 | 0/1 | 1.20.1440.150 2.40.50.140 |
| 2jjlA | 2.30 | 207 | 944 | 3 | 35 40 132 | 1.000 | 9.12 | 0.582 | 805 | 7.78 | 0.574 | 0/4 | 1.20.5.170 2.10.25.20 2.60.90.40 |
| 1ro2A | 1.60 | 210 | 1003 | 3 | 70 115 25 | 1.000 | 9.55 | 0.569 | 840 | 8.00 | 0.573 | 0/6 | 2.20.25.130 3.30.2250.10 3.30.720.160 |
| 2x49A | 1.50 | 333 | 1562 | 4 | 96 48 90 99 | 1.000 | 9.38 | 0.559 | 1295 | 7.78 | 0.580 | 0/0 | 1.10.8.540 3.40.30.60 3.40.5.40 3.40.50.12790 |
| 1h3jB | 2.00 | 336 | 1702 | 2 | 162 162 | 0.964 | 10.13 | 0.518 | 1323 | 7.88 | 0.523 | 0/6 | 1.10.420.10 1.10.520.10 |
| 7k76L | 2.14 | 213 | 1036 | 2 | 102 97 | 0.934 | 9.73 | 0.540 | 836 | 7.85 | 0.538 | 0/0 | 2.60.40.10 2.60.40.10 |
| 2a7kA | 2.24 | 229 | 1186 | 2 | 196 33 | 1.000 | 10.36 | 0.544 | 968 | 8.45 | 0.566 | 1/21 | 1.20.5.1610 3.90.226.10 |
| 3ihmB | 2.30 | 407 | 2063 | 3 | 214 144 49 | 1.000 | 10.14 | 0.532 | 1659 | 8.15 | 0.537 | 0/23 | 3.30.9.40 3.50.50.60 6.10.250.650 |
| 1rysB | 2.03 | 341 | 1678 | 4 | 67 98 67 93 | 0.953 | 9.84 | 0.559 | 1370 | 8.04 | 0.566 | 0/11 | 1.10.150.20 3.30.1490.100 3.30.70.270 3.40.1170.60 |
| 3cw7D | 2.30 | 282 | 1351 | 3 | 112 118 52 | 1.000 | 9.58 | 0.546 | 1102 | 7.82 | 0.559 | 0/0 | 1.10.1670.10 1.10.340.30 3.30.310.20 |
| 1ulyA | 2.50 | 190 | 863 | 2 | 89 101 | 1.000 | 9.08 | 0.578 | 739 | 7.78 | 0.594 | 0/2 | 1.10.10.10 1.10.20.10 |
| 6v06A | 2.40 | 326 | 1540 | 5 | 62 59 62 60 82 | 0.997 | 9.45 | 0.544 | 1218 | 7.47 | 0.535 | 0/0 | 2.10.70.10 2.10.70.10 2.10.70.10 2.10.70.10 2.10.70.10 |
| 1gz0A | 2.50 | 242 | 1164 | 2 | 77 165 | 1.000 | 9.62 | 0.556 | 944 | 7.80 | 0.558 | 0/11 | 3.30.1330.30 3.40.1280.10 |
| 1b7yB | 2.50 | 775 | 3731 | 6 | 70 116 214 74 191 86 | 0.969 | 9.63 | 0.543 | 3016 | 7.78 | 0.546 | 0/10 | 2.40.50.140 3.30.56.10 3.30.56.10 3.30.70.380 3.30.930.10 3.50.40.10 |
| 3fulA | 2.39 | 605 | 3085 | 4 | 193 97 157 148 | 0.983 | 10.20 | 0.522 | 2462 | 8.14 | 0.539 | 0/6 | 1.10.390.10 1.25.40.320 2.60.40.1730 3.30.2010.30 |
| 6u45A | 2.35 | 726 | 3630 | 4 | 249 112 76 153 | 0.813 | 10.00 | 0.538 | 2928 | 8.07 | 0.549 | 0/6 | 2.40.30.10 3.30.230.10 3.30.70.870 3.40.50.300 |
| 4ihfA | 2.10 | 336 | 1798 | 2 | 97 239 | 1.000 | 10.70 | 0.557 | 1450 | 8.63 | 0.532 | 0/12 | 2.160.10.10 3.40.1390.10 |
| 3sluB | 2.41 | 358 | 1769 | 3 | 91 138 129 | 1.000 | 9.88 | 0.548 | 1422 | 7.94 | 0.547 | 5/13 | 2.70.70.10 3.10.450.350 3.10.450.350 |
| 2bnoA | 1.90 | 191 | 901 | 2 | 75 116 | 1.000 | 9.43 | 0.581 | 751 | 7.86 | 0.567 | 4/7 | 1.10.260.40 2.60.120.10 |
| 4c73B | 1.80 | 414 | 2355 | 2 | 264 150 | 1.000 | 11.38 | 0.522 | 1823 | 8.81 | 0.517 | 0/25 | 3.40.47.10 3.40.47.10 |
| 1qupA | 1.80 | 221 | 1068 | 2 | 70 151 | 1.000 | 9.67 | 0.547 | 861 | 7.79 | 0.549 | 0/1 | 2.60.40.200 3.30.70.100 |
| 2d1cA | 1.80 | 495 | 2520 | 2 | 384 111 | 1.000 | 10.18 | 0.534 | 1999 | 8.08 | 0.537 | 0/1 | 3.30.70.1570 3.40.718.10 |
| 2zo9B | 2.20 | 271 | 1356 | 2 | 131 126 | 0.948 | 10.01 | 0.537 | 1074 | 7.93 | 0.532 | 0/3 | 3.30.750.180 3.60.21.40 |
| 1nxzB | 2.00 | 245 | 1149 | 2 | 72 173 | 1.000 | 9.38 | 0.545 | 958 | 7.82 | 0.541 | 0/1 | 2.40.240.20 3.40.1280.10 |
| 1cyyA | 2.15 | 250 | 1188 | 2 | 125 125 | 1.000 | 9.50 | 0.555 | 987 | 7.90 | 0.556 | 5/14 | 1.10.290.10 2.70.20.10 |
| 1v6vB | 2.10 | 436 | 2297 | 2 | 122 300 | 0.968 | 10.54 | 0.529 | 1816 | 8.33 | 0.525 | 0/0 | 2.80.10.50 3.20.20.80 |
| 3smzA | 1.99 | 280 | 1359 | 3 | 112 71 97 | 1.000 | 9.71 | 0.555 | 1106 | 7.90 | 0.557 | 0/4 | 3.30.70.330 3.30.70.330 3.30.70.330 |
| 1xdpA | 2.50 | 687 | 3315 | 4 | 105 202 188 185 | 0.990 | 9.65 | 0.543 | 2736 | 7.97 | 0.551 | 0/0 | 1.20.58.310 3.30.1840.10 3.30.870.10 3.30.870.10 |
| 3mogA | 2.20 | 470 | 2399 | 3 | 187 189 94 | 1.000 | 10.21 | 0.532 | 1958 | 8.33 | 0.539 | 0/13 | 1.10.1040.50 3.30.750.190 3.40.50.720 |
| 5le8B | 1.78 | 318 | 1628 | 2 | 154 163 | 0.997 | 10.24 | 0.535 | 1342 | 8.44 | 0.536 | 0/12 | 1.25.40.20 1.25.40.20 |
| 1x0cB | 1.70 | 549 | 2972 | 2 | 171 374 | 0.993 | 10.83 | 0.516 | 2295 | 8.36 | 0.502 | 0/0 | 2.160.20.10 2.60.350.10 |
| 3m7nB | 2.40 | 179 | 880 | 3 | 52 86 34 | 0.961 | 9.83 | 0.561 | 696 | 7.78 | 0.553 | 0/0 | 2.20.70.10 2.40.50.100 2.40.50.140 |
| 2nv9F | 1.95 | 369 | 1861 | 2 | 127 231 | 0.970 | 10.09 | 0.537 | 1505 | 8.16 | 0.540 | 3/3 | 2.40.37.10 3.20.20.10 |
| 5z24A | 2.40 | 441 | 2173 | 3 | 145 175 118 | 0.993 | 9.85 | 0.547 | 1746 | 7.92 | 0.535 | 3/24 | 2.60.40.10 2.60.40.10 2.60.40.740 |
| 1it8A | 2.50 | 577 | 3010 | 4 | 355 70 74 78 | 1.000 | 10.43 | 0.519 | 2389 | 8.28 | 0.539 | 0/5 | 2.30.130.10 3.10.450.90 3.20.20.105 3.90.1020.10 |
| 4yb5D | 2.24 | 292 | 1385 | 3 | 131 87 74 | 1.000 | 9.49 | 0.570 | 1131 | 7.75 | 0.569 | 3/8 | 3.30.70.120 3.40.190.10 3.40.190.10 |
| 2fjrB | 1.95 | 189 | 841 | 2 | 73 116 | 1.000 | 8.90 | 0.580 | 700 | 7.41 | 0.578 | 0/0 | 1.10.260.40 2.10.109.10 |
| **mean over 60 chains** | |369.4 | 1852.9 | 2.78 | | 0.983 | 9.92 | 0.545 | 1497.7 | 8.05 | 0.549 | | |

## Caveats

- CATH chops domains partly from 3D compactness and contacts, so domains are
  expected to be dense in a contact graph; the task is easier than ground truth
  that is independent of the edges.
- The 8 A C-alpha rule joins every residue to its chain neighbours, so the
  graphs carry a strong 1D backbone signal; most domains are one or two
  contiguous sequence segments.
- Only one chain of each crystal is used; contacts with other chains,
  ligands or crystal mates are ignored.
- The missing-residue and coverage filters favour well-ordered, fully
  classified chains.
- Sampling is uniform over superfamily combinations, not over chains or
  over the PDB.
