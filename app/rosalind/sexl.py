


from typing import List


def female(q: float):
    p = 1 - q
    return 2 * p * q

def sex_inheritance(males: List[float]) -> List[float]:
    females = [female(i) for i in males]
    return females
