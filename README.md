# HyperCommon

HyperCommon finds overlapping communities in networks. Three nodes, one adjacent to the
other two, whose three pairs all reach a commonality (neighbourhood-overlap) threshold
form an elementary community; communities are the connected components of the resulting
3-uniform hypergraph under shared pairs. The construction is deterministic, local and
nested in its threshold.

This repository holds the code and the full reproduction pipeline of the paper

> D. Deinega, I. Kanovsky and B. Barzel, *Overlapping communities as connected components
> of a commonality hypergraph* (submitted to Physical Review E).

**To reproduce the paper, start at [`paper/README.md`](paper/README.md).** It lists the
software versions, the data sources, the seeds, the hardware of every timing and the
command behind every table and figure.

## Layout

| Path | What it is |
|---|---|
| `hypercommon/` | reference implementation (hypergraph construction, communities, local query) |
| `predicates/`, `metrics/`, `generators/` | Jaccard commonality, the Omega index, ring-lattice generators |
| `tests/` | tests of the reference implementation (`python -m pytest tests`) |
| `paper/` | the manuscript, the benchmark and analysis pipeline, the raw results and the vectorised implementation used for the timings |
| `experiments/ring_lattice/`, `plots/`, `results/`, `main.py` | earlier exploratory ring-lattice experiments; the paper does not use them |
| `experiments/gnn/`, `results/gnn/` | an unrelated exploratory experiment; the paper does not use it |

## Install

```bash
python -m pip install -r requirements.txt
```

Python 3.12. The stochastic block model baseline needs graph-tool under Linux
(`paper/environment-sbm.yml`).

## Licence

MIT (see `LICENSE`). If you use the code, please cite the paper (`CITATION.cff`).
