

from typing import Tuple


def transition_transversion_count(seq1: str, seq2: str) -> Tuple[int, int]:
    itions = 0
    versions = 0
    for a, b in zip(seq1, seq2):
        if not a == b:
            if a in "AG" and b in "AG":
                itions += 1
            elif a in "CT" and b in "CT":
                itions += 1
            else:
                versions += 1
    return itions, versions

