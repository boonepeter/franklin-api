

from typing import List, Set


def set_operations(size: int, set1: set, set2: set) -> List[Set[int]]:
    big = set(i + 1 for i in range(size))
    union = set1.union(set2)
    inter = set1.intersection(set2)
    not1 = set1.difference(set2)
    not2 = set2.difference(set1)
    comp1 = big.difference(set1)
    comp2 = big.difference(set2)

    return [union, inter, not1, not2, comp1, comp2]