from fastapi import HTTPException

COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

def reverse_complement(dna: str) -> str:
    complement = []
    for b in dna[::-1]:
        if b not in COMPLEMENT:
            raise HTTPException(status_code=500, detail=f"{b} not in {COMPLEMENT.keys()}")
        complement.append(COMPLEMENT[b])
    return "".join(complement)