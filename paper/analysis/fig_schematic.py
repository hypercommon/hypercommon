"""
fig_schematic.pdf: the method on a toy graph, computed by the actual code.

(a) Two ring-lattice groups (k=4) sharing one merged node, plus one joining edge
    between them. Admitted triples at threshold t are drawn as translucent
    triangles coloured by their component; the shared node belongs to both
    communities; the joining edge (4-14, closing no triangle with the shared node)
    carries no admitted triple.
(b) The two-node-overlap graph L_t of the admitted triples: one vertex per
    triple, an edge when two triples share a pair. Its components are the
    communities.
"""
import os
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench"))
from hc_fast import TripleScores                     # noqa: E402
from fig_real_networks import INK, INK2, FIG          # noqa: E402

T = 0.25
COLS = ["#2a78d6", "#eb6834"]

# two rings of 9, each node linked to 2 neighbours per side; node 0 of ring B merged into node 0 of ring A
G = nx.Graph()
for r, off in ((0, 0), (1, 9)):
    for i in range(9):
        for d in (1, 2):
            G.add_edge(off + i, off + (i + d) % 9)
for nb in list(G.neighbors(9)):
    G.add_edge(0, nb)
G.remove_node(9)
G.add_edge(4, 14)            # a joining edge between the groups, far from the shared node

ts = TripleScores(G)
comms = ts.communities("hc", T)
comm_of_ring_a = next(i for i, C in enumerate(comms) if 3 in C)
keep = np.flatnonzero(ts.s3 >= T)
tri = ts.triples[keep]
# component of each admitted triple = the community containing all three nodes
nodes = ts.nodes
tri_comm = []
for a, b, c in tri:
    s = {nodes[a], nodes[b], nodes[c]}
    tri_comm.append(next(i for i, C in enumerate(comms) if s <= C))

# layout: two circles side by side, shared node between them
pos = {0: (0.0, 0.0)}                      # the shared node, where the rings touch
for i in range(1, 9):                       # ring A around (-1, 0), node i at angle 2*pi*i/9
    ang = 2 * np.pi * i / 9
    pos[i] = (-1 + np.cos(ang), np.sin(ang))
for i in range(1, 9):                       # ring B around (+1, 0), starting from the left
    ang = np.pi + 2 * np.pi * i / 9
    pos[9 + i] = (1 + np.cos(ang), np.sin(ang))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.35, 1]})
for (a, b, c), ci in zip(tri, tri_comm):
    pts = [pos[nodes[a]], pos[nodes[b]], pos[nodes[c]]]
    ax.add_patch(Polygon(pts, closed=True, color=COLS[ci], alpha=0.07, lw=0))
for u, v in G.edges():
    random_edge = {u, v} == {4, 14}      # the joining edge
    ax.plot(*zip(pos[u], pos[v]), color=INK2 if not random_edge else "#e34948",
            lw=0.8 if not random_edge else 1.2, ls="-" if not random_edge else (0, (3, 2)), zorder=1)
for u in G.nodes():
    member = [i for i, C in enumerate(comms) if u in C]
    if len(member) == 2:
        ax.add_patch(plt.matplotlib.patches.Wedge(pos[u], 0.13, 90, 270, color=COLS[0], zorder=3))
        ax.add_patch(plt.matplotlib.patches.Wedge(pos[u], 0.13, -90, 90, color=COLS[1], zorder=3))
    else:
        ax.add_patch(plt.Circle(pos[u], 0.11, color=COLS[member[0]] if member else "white",
                                ec=INK, lw=0.6, zorder=3))
ax.annotate("shared node\n(both communities)", xy=(0, -0.13), xytext=(0, -1.05), ha="center", va="center",
            fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.6))
ax.text(0, 1.12, "joining edge: no admissible triple", ha="center", fontsize=6.3, color="#e34948")
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.2, 1.3); ax.set_aspect("equal"); ax.axis("off")
ax.set_title(rf"(a) admitted triples at $t={T}$", fontsize=8.5, loc="left", fontweight="bold")

# (b) the shared-pair graph of admitted triples
L = nx.Graph()
pair_owner = {}
for ti, (a, b, c) in enumerate(tri):
    L.add_node(ti)
    for p in ((a, b), (a, c), (b, c)):
        for tj in pair_owner.get(p, []):
            L.add_edge(ti, tj)
        pair_owner.setdefault(p, []).append(ti)
# each component on its own circle, a triple placed at the angle of its centroid around
# its ring's centre in panel (a); triples through the shared node are ringed in black
lp = {}
for ti, (a, b, c) in enumerate(tri):
    ci = tri_comm[ti]
    cx = -1.0 if ci == comm_of_ring_a else 1.0
    px = np.mean([pos[nodes[x]][0] for x in (a, b, c)]) - cx
    py = np.mean([pos[nodes[x]][1] for x in (a, b, c)])
    ang = np.arctan2(py, px)
    rad = 1.0 + 0.12 * ((ti * 7) % 3 - 1)
    lp[ti] = (1.25 * cx + rad * np.cos(ang), rad * np.sin(ang))
for u, v in L.edges():
    bx.plot(*zip(lp[u], lp[v]), color=INK2, lw=0.35, alpha=0.55, zorder=1)
for ti, (a, b, c) in enumerate(tri):
    through = 0 in (nodes[a], nodes[b], nodes[c])
    bx.scatter(*lp[ti], s=22 if through else 14, color=COLS[tri_comm[ti]], zorder=2,
               edgecolors=INK if through else "white", linewidths=0.9 if through else 0.4)
bx.text(0, -1.35, "black rim: triples through the shared node", ha="center", fontsize=6.3, color=INK2)
bx.set_xlim(-2.6, 2.6); bx.set_ylim(-1.45, 1.3); bx.set_aspect("equal")
bx.axis("off")
bx.set_title(r"(b) the graph $L_t$ of admitted triples", fontsize=8.5, loc="left", fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig_schematic.pdf"))
fig.savefig(os.path.join(FIG, "fig_schematic.png"), dpi=220)
print("communities:", [sorted(c) for c in comms], "admitted triples:", len(tri))
