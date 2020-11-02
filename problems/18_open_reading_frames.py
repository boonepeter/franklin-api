"""
Problem
Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

Given: A DNA string s of length at most 1 kbp in FASTA format.

Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.

Sample Dataset
>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
Sample Output
MLLGSFRLIPKETLIQVAGSSPCNLS
M
MGMTPRLGLESLLE
MTPRLGLESLLE
"""
from typing import List
from utils.read_fasta import read_fasta
from utils.codons import rna_to_protein

filename = "./input/rosalind_orf.txt"
seq = read_fasta(filename)

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
        



for s in seq.values():
    rna = s.replace("T", "U")
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
    for o in set(orfs):
        print(o)