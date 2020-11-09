

from typing import List




"""
d(seq) = min number of reversals
b(seq) = number of breakpoints

d(seq) >= b(seq) / 2

reversal is sorting if it reduces the distance to identity by 1

reversal can remove 2 breakpoints without being sorting
"""

def num_breakpoints(seq: List[int]) -> int:
    bp_count = 0
    n_list = [min(seq) - 1] + seq + [max(seq) + 1]
    for i in range(1, len(n_list)):
        if not abs(n_list[i] - n_list[i - 1]) == 1:
            bp_count += 1
    return bp_count

def sort(seq: List[int]) -> int:

    return 0


def transform(seq1: List[int], seq2: List[int], offest: int=1) -> List[int]:
    """
    Orders seq2 by seq1. Makes this just a sorting problem
    """
    t = [0 for i in seq1]
    for i, n in enumerate(seq1):
        t[n - offest] = seq2[i]
    return t


    