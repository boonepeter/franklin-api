from fastapi import HTTPException
from ..utils.codons import protein_to_rna

def infer_rna_combinations(protein: str, mod: int=1000000) -> int:
    total = 1
    for p in protein:
        if p not in protein_to_rna:
            raise HTTPException(status_code=500, detail=f"{p} not in codon table")
        total *= len(protein_to_rna[p])
        total %= mod
    total *= 3
    total %= mod
    return total