
from math import factorial
def perfect_rna_match(rna: str) -> int:
    return factorial(rna.count('A')) * factorial(rna.count('C'))