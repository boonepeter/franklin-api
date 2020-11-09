
from itertools import permutations
from typing import List, Tuple

def signed_probs(num: int) -> List[Tuple[int, ...]]:
    nums = [(i, -i) for i in range(1, num + 1)]
    nums = [i for sublist in nums for i in sublist]
    perms = list(permutations(nums, r=num))
    good = []
    for p in perms:
        bad = False
        for i in range(1, num + 1):
            if i in p and -i in p:
                bad = True
                break
        if not bad:
            good.append(p)
    return good