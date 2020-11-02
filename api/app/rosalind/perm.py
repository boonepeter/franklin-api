
from itertools import permutations
from typing import List, Tuple

def all_permutations(num: int) -> List[Tuple[int, ...]]:
    """
    Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).
    """
    perms = [i for i in permutations(range(1, num + 1))]
    return perms