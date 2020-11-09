
from typing import List, Tuple


REV = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

def is_palindrom(seq: str) -> bool:
    for i in range(len(seq)):
        a = seq[i]
        b = seq[-(i + 1)]
        if not a == REV[b]:
            return False
    return True

def get_palindrom_locations(sequence: str, min_len: int=4, max_len: int=12) -> List[Tuple[int, int]]:
    pals = []
    for i in range(len(sequence)):
        for j in range(min_len, max_len + 1):
            if i + j > len(sequence):
                continue
            s = sequence[i: i + j]
            if is_palindrom(s):
                pals.append((i + 1, len(s)))
    return pals