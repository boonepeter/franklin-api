
from math import factorial

def combination_n_k(n: int, k: int) -> int:
    return factorial(n) // (factorial(k) * factorial(n - k))


def alt_splicing(n: int, m: int, mod:int=1000000) -> int:
    total = 0
    for k in range(m, n + 1):
        total += combination_n_k(n, k)
        total = total % mod
    return total