from typing import List


def find_motifs(sequence: str, motif: str) -> List[int]:
    positions = []
    for i in range(0, len(sequence) - len(motif)):
        if sequence[i: i + len(motif)] == motif:
            positions.append(i + 1)
    return positions