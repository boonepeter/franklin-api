"""
Problem
After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.

Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)

Sample Dataset
>Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT
Sample Output
MVYIADKQHVASREAYGHMFKVCA
"""
from utils.read_fasta import read_fasta
from utils.codons import translate
filename = "./input/rosalind_splc.txt"

seqs = read_fasta(filename)

items = list(seqs.values())

DNA: str = items[0]

for i in items[1:]:
    DNA = DNA.replace(i, "")

print(translate(DNA))
