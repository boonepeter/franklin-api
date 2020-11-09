

from typing import List

def distance(seq1: str, seq2: str) -> float:
    dont_match = 0
    for a, b in zip(seq1, seq2):
        if not a == b:
            dont_match += 1
    return dont_match / len(seq1)


def p_distance(seqs: List[str]) -> List[List[float]]:
    dists = [[0.0 for i in seqs] for i in seqs]

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i == j:
                continue
            if dists[i][j] > 0:
                dists[j][j] = dists[i][j]
                continue
            elif dists[j][i] > 0:
                dists[i][j] = dists[j][i]
            else:
                dists[i][j] = distance(seqs[i], seqs[j])
    return dists