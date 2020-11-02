from typing import Dict
from fastapi import HTTPException

def count_dna(seq: str) -> Dict[str, int]:
    counts = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0
    }
    for b in seq.upper():
        if b not in counts:
            raise HTTPException(status_code=500, detail=f"{b} not in {counts.keys()}")
        counts[b] += 1
    return counts