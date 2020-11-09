from typing import List

def remove_introns(first: str, introns: List[str]) -> str:
    for i in introns:
        first = first.replace(i, "")
    return first