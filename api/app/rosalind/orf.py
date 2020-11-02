
from typing import List, Set
from ..utils.codons import rna_to_protein

def all_orf(seq: str, start) -> List[str]:
    potential = []
    prev = ""
    running: List[str] = []
    for i in range(start, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) < 3:
            break
        protein = rna_to_protein[codon]
        if protein == "Stop":
            if len(running) > 0:
                potential.extend(running)
            running = []
            continue
        elif protein == "M":
            running = [r + protein for r in running]
            running.append(protein)
            continue
        if len(running) > 0:
            running = [r + protein for r in running]    
    return potential
        

def get_all_orfs(seq: str) -> Set[str]:
    rna = seq.replace("T", "U")
    orfs = []

    complement = ""
    comp = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C"
    }

    for b in rna[::-1]:
        complement += comp[b]

    orfs.extend(all_orf(rna, 0))
    orfs.extend(all_orf(rna, 1))
    orfs.extend(all_orf(rna, 2))
    orfs.extend(all_orf(complement, 0))
    orfs.extend(all_orf(complement, 1))
    orfs.extend(all_orf(complement, 2))
    return set(orfs)

