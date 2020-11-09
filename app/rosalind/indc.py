

from typing import List

from math import factorial, log

def binomial(n: int, k: int, p: float) -> float:
    f = factorial(n) / factorial(n - k) / factorial(k)
    return f * p**k * (1 - p) ** (n - k)

def ind_segregation(n: int, p: float=0.5) -> List[float]:
    vals = []
    for k in range(1, 2 * n + 1):
        bins = [binomial(2 * n, kk, p) for kk in range(k, 2 * n + 1)]
        vals.append(log(sum(bins), 10))
    return vals