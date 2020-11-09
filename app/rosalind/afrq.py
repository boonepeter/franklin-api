

from typing import List


def recessive_probability(homo_pop: List[float], p: float=0.5) -> List[float]:
    """homo_pop is the portion of the pop that are homozygous.
    We want how many will have at least one allele.
    """
    probs: List[float] = []
    for i in homo_pop:
        prob = -i + 2 * i ** p
        probs.append(prob)
    return probs