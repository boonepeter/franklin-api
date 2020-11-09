

from app.rosalind.indc import binomial
from scipy.stats.distributions import binom

def wright_fisher(N: int, m: int, g: int, k: int) -> float:
    p = m / (2 * N)
    q = 1 - p
    bc = binomial(N, k, p)
    print(binom.cdf(k, N, p))
    return bc