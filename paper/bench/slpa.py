"""SLPA, as methods2._slpa_multi (bundled for big_job.py)."""
import numpy as np


def slpa_multi(G, T, rs):
    """SLPA (Xie, Szymanski & Liu 2011) as in cdlib's slpa_nx, with the label-propagation
    stage run once and its post-processing threshold r applied for every r of the grid
    (the threshold only filters the final memories)."""
    from collections import defaultdict
    memory = {i: {i: 1} for i in G.nodes()}
    for _ in range(T):
        order = list(G.nodes())
        np.random.shuffle(order)
        for listener in order:
            speakers = G[listener].keys()
            if len(speakers) == 0:
                continue
            labels = defaultdict(int)
            for speaker in speakers:
                mem = memory[speaker]
                total = float(sum(mem.values()))
                labels[list(mem.keys())[np.random.multinomial(1, [f / total for f in mem.values()]).argmax()]] += 1
            acc = max(labels, key=labels.get)
            memory[listener][acc] = memory[listener].get(acc, 0) + 1
    for r in rs:
        comms = {}
        for node, mem in memory.items():
            for label, freq in mem.items():
                if freq / float(T + 1) >= r:
                    comms.setdefault(label, set()).add(node)
        keys = list(comms)
        nested = set()
        for i, a in enumerate(keys[:-1]):
            for b in keys[i + 1:]:
                if comms[a] <= comms[b]:
                    nested.add(a)
                elif comms[a] >= comms[b]:
                    nested.add(b)
        yield r, [c for k, c in comms.items() if k not in nested and len(c) > 0]
