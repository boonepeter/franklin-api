

from scipy.special import binom


def independent_alleles(k: int, N: int) -> int:
    def p(n: int, k: int):
        return binom(2 ** k, n) * 0.25 ** n * 0.75 ** (2 ** k - n)
    return 1 - sum(p(n, k) for n in range(N))
