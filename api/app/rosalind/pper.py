
from math import factorial
def partial_permutations(n: int, k: int, mod: int=1000000) -> int:
    parts = (factorial(n) // factorial(n - k)) % mod
    return parts