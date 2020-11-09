
_table = """UUU F      CUU L      AUU I      GUU V
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
rna_to_protein = {}
for l in _table.split("\n"):
    splits = l.split()
    for i in range(0, len(splits), 2):
        rna_to_protein[splits[i]] = splits[i + 1]


protein_to_rna = {}
for k, v in rna_to_protein.items():
    if v not in protein_to_rna:
        protein_to_rna[v] = [k]
    else:
        protein_to_rna[v].append(k)

def translate(seq: str):
    seq = seq.replace("T", "U")
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        aa = rna_to_protein[seq[i: i + 3]]
        if aa == "Stop":
            return protein
        protein += aa
    return protein
