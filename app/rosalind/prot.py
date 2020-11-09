from fastapi import HTTPException

table = """UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G"""
codons = {}
for l in table.split("\n"):
    splits = l.split()
    for i in range(0, len(splits), 2):
        codons[splits[i]] = splits[i + 1]


def rna_to_protein(rna: str, with_stop=False) -> str:
    protein = ""
    for i in range(0, len(rna), 3):
        codon = rna[i: i + 3].upper()
        if len(codon) != 3:
            continue
        if codon not in codons:
            raise HTTPException(status_code=500, detail=f"{codon} not in the codon table")
        aa = codons[rna[i: i + 3]]
        if aa == "Stop" and not with_stop:
            continue
        protein += aa
    return protein

def dna_to_protein(dna: str, with_stop=False) -> str:
    rna = dna.replace("T", "U")
    return rna_to_protein(rna, with_stop)
