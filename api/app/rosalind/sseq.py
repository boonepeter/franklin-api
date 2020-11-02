

from typing import List

def splice_motifs(dna: str, motif: str) -> List[int]:
    indices = []
    for i, bp in enumerate(dna):
        if len(motif) == 0:
            break
        if motif[0] == bp:
            indices.append(i + 1)
            motif = motif[1:]
    return indices