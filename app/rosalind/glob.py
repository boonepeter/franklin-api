

from ..utils.blosum62 import cost as blosum_cost
from itertools import product

def glob_alignment(s: str, t: str, gap_penalty:int=5) -> int:
    sl, tl = len(s), len(t)
    m = {(0, 0): (0, None)}
    m.update({((i, 0), (i * -5, (i-1, 0))) for i in range(1, sl+1)})
    m.update({((0, i), (i * -5, (0, i-1))) for i in range(1, tl+1)})
    i, j = 0, 0
    for i, j in product(range(1, sl+1), range(1, tl+1)):
        cost = blosum_cost(s[i-1], t[j - 1])
        d = m[(i-1, j-1)][0] + cost
        l = m[(i-1, j)][0] - gap_penalty
        u = m[(i, j-1)][0] - gap_penalty
        b = max(d, l, u)
        if d == b:
            m[(i, j)] = (b, (i - 1, j - 1))
        elif l == b:
            m[(i, j)] = (b, (i - 1, j))
        elif u == b:
            m[(i, j)] = (b, (i, j - 1))
    return m[(i, j)][0]