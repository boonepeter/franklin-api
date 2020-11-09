
from typing import List
from math import log10

def probability(seq: str, nums: List[float]) -> List[float]:

    bp_count = {
        "G": 0,
        "C": 0,
        "A": 0,
        "T": 0
    }

    for bp in seq:
        bp_count[bp] += 1

    gc = bp_count["G"] + bp_count["C"]
    at = bp_count["A"] + bp_count["T"]

    gc_probs = []
    for v in nums:
        log_prob = gc * log10(0.5 * v) + at * log10(0.5 * (1 - v))
        gc_probs.append(log_prob)
    
    return gc_probs