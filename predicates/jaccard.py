from hypercommon.hypernode import HCNode


def closed_neighborhood_jaccard(u: HCNode, v: HCNode) -> float:
    Nu = u.neighbors | {u.id}
    Nv = v.neighbors | {v.id}
    return len(Nu & Nv) / len(Nu | Nv)


def closed_neighborhood_jaccard_predicate(threshold: float):
    def predicate(u: HCNode, v: HCNode) -> bool:
        return closed_neighborhood_jaccard(u, v) >= threshold
    return predicate