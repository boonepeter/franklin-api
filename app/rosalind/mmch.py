
from math import factorial

def permutations(n: int, k: int):
    return factorial(n) // factorial(n-k)


def num_matches(seq: str):
    na = seq.count('A')
    nu = seq.count('U')
    nc = seq.count('C')
    ng = seq.count('G')
    max_au = max(na, nu)
    min_au = min(na, nu)

    max_gc = max(nc, ng)
    min_gc = min(nc, ng)


    return (permutations(max_au, min_au) *
            permutations(max_gc, min_gc))