from typing import List, Tuple
from itertools import product

def lex_product(letters: List[str], length: int) -> List[Tuple[str, ...]]:
    letters.sort()
    prods = product(letters, repeat=length)
    return [i for i in prods]

