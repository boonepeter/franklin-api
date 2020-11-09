
from typing import List
from itertools import product

def lex_permutations(alphabet: List[str], length: int, to_sort: bool=True):
    prods = []
    for i in range(length):
        prods.extend([p for p in product(alphabet, repeat=i + 1)])
    joined = ["".join(p) for p in prods]
    if to_sort:
        joined = sorted(joined, key=lambda word: [alphabet.index(c) for c in word])
    return joined
