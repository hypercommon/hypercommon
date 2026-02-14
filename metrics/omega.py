from collections import defaultdict
from itertools import combinations


def build_pair_counts(communities):
    pair_counts = defaultdict(int)
    for community in communities:
        for i, j in combinations(community, 2):
            if i > j:
                i, j = j, i
            pair_counts[(i, j)] += 1
    return pair_counts


def omega_index(
    ground_truth_pair_counts,
    detected_pair_counts,
    total_pairs: int,
):
    """
    Compute the adjusted Omega Index for overlapping community detection.

    Input partitions are represented by pair co-membership counts:
      - ground_truth_pair_counts[(i,j)] = how many ground-truth communities contain both i and j
      - detected_pair_counts[(i,j)]     = how many detected communities contain both i and j

    Pairs not present in a dict are assumed to have co-membership count 0.

    The Omega Index measures agreement between the two partitions over all unordered node pairs,
    and adjusts for the agreement expected by chance:
      - Omega_u: fraction of node pairs whose co-membership counts match exactly
      - Omega_e: expected Omega_u under independence (computed from histograms of counts)
      - Omega   = (Omega_u - Omega_e) / (1 - Omega_e)

    Parameters
    ----------
    ground_truth_pair_counts : Mapping[(int,int), int]
        Co-membership counts for the ground truth partition.
    detected_pair_counts : Mapping[(int,int), int]
        Co-membership counts for the detected partition.
    total_pairs : int
        Total number of unordered node pairs in the evaluated node universe (N*(N-1)//2).

    Returns
    -------
    float
        Adjusted Omega Index in [-1, 1] (1 = perfect agreement, 0 = chance-level agreement).
    """

    freq_gt = defaultdict(int)
    freq_det = defaultdict(int)

    # pairs that are nonzero in at least one partition
    keys = set(ground_truth_pair_counts.keys()) | set(detected_pair_counts.keys())

    agree_nonzero = 0

    for key in keys:
        c1 = ground_truth_pair_counts.get(key, 0)
        c2 = detected_pair_counts.get(key, 0)

        if c1 == c2:
            agree_nonzero += 1

        freq_gt[c1] += 1
        freq_det[c2] += 1

    # remaining pairs are (0,0)
    remaining = total_pairs - len(keys)
    if remaining < 0:
        raise ValueError("total_pairs is smaller than number of observed keys. Wrong N/total_pairs?")

    # all those remaining pairs agree (both 0)
    agree_total = agree_nonzero + remaining

    # add their histogram mass
    freq_gt[0] += remaining
    freq_det[0] += remaining

    Omega_u = agree_total / total_pairs

    Omega_e = 0.0
    for k in set(freq_gt.keys()) | set(freq_det.keys()):
        p1 = freq_gt[k] / total_pairs
        p2 = freq_det[k] / total_pairs
        Omega_e += p1 * p2

    if Omega_e == 1.0:
        return 1.0
    return (Omega_u - Omega_e) / (1.0 - Omega_e)