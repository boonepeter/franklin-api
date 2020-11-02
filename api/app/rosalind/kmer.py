from typing import DefaultDict, List, Dict
from itertools import product

def kmer_count(seq: str, k:int) -> Dict[str, int]:
    counts = {}
    a = ["".join(i) for i in product("ACGT", repeat=k)]

    for i in a:
        counts[i] = 0
    for i in range(0, len(seq) - k + 1):
        sub = seq[i: i + 4]
        counts[sub] += 1
    return counts

def kmer_array(seq: str, k: int=4) -> List[int]:
    counts = kmer_count(seq, k)
    keys = list(counts.keys())
    keys.sort()
    arr: List[int] = []
    for key in keys:
        arr.append(counts[key])
    return arr
    

