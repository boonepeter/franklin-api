

def hamming_distance(one: str, two: str) -> int:
    distance = 0
    for o, t in zip(one, two):
        if not o == t:
            distance += 1
    return distance