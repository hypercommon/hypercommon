import networkx as nx
import pandas as pd

from generators.ring_lattice import ring_lattice
from hypercommon.hypergraph import build_hypergraph
from predicates.jaccard import closed_neighborhood_jaccard_predicate
from utils.threshold import representative_thresholds


def degree_stats(H: nx.Graph):
    if H.number_of_nodes() == 0:
        return {
            "Degree Set": "{}"
        }

    degs = [d for _, d in H.degree]
    s = sorted(set(degs))

    return {
        "Degree Set": "{" + ",".join(map(str, s)) + "}"
    }


# --------------------------------------------------
# Experiment runner
# --------------------------------------------------

# configs: iterable of (n, z) pairs (rings is fixed to 1 in this experiment)
def run_experiment(configs, output_path="results/commonality_thresholds.xlsx"):
    writer = pd.ExcelWriter(output_path, engine="openpyxl")

    for n, z in configs:
        rings = 1

        # Generate lattice
        G = ring_lattice(n=n, z=z, rings=rings)
        E = G.number_of_edges()

        # Thresholds
        thresholds = representative_thresholds(G)

        rows = []

        for t in thresholds:
            H = build_hypergraph(G, closed_neighborhood_jaccard_predicate(t))

            stats = degree_stats(H)

            rows.append({
                "Threshold": round(t, 6),
                "H Nodes": H.number_of_nodes(),
                "H Edges": H.number_of_edges(),
                **stats
            })

        df = pd.DataFrame(rows)

        sheet_name = f"n{n}_e{E}_z{z}"
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    writer.close()