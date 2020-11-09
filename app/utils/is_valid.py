
from fastapi import HTTPException



def raise_valid_dna(seq: str):
    """Raise an HTTPException if an invalid character is found (not in ACGT)"""
    if not is_valid_dna(seq):
        raise HTTPException(status_code=500, detail="Sequence is not valid DNA")

def raise_valid_rna(seq: str):
    """Raise an HTTPException is an invalid character is found (not in ACGU)"""
    if not is_valid_rna(seq):
        raise HTTPException(status_code=500, detail="Not a valid RNA sequence")

def is_valid_dna(seq: str) -> bool:
    for bp in seq:
        if bp not in "ACGT":
            return False
    return True

def is_valid_rna(seq: str) -> bool:
    for bp in seq:
        if bp not in "ACGU":
            return False
    return True
