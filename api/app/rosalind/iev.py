
from typing import List


PROB = [1, 1, 1, 0.75, 0.5, 0]

def expected(counts: List[int], number=2) -> float:
    total = 0
    for p, c in zip(PROB, counts):
        total += c * p * number
    return total