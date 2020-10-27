"""
Problem
To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.

You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into

http://www.uniprot.org/uniprot/uniprot_id
Alternatively, you can obtain a protein sequence in FASTA format by following

http://www.uniprot.org/uniprot/uniprot_id.fasta
For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.

Given: At most 15 UniProt Protein Database access IDs.

Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.

Sample Dataset
A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST
Sample Output
B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614
Noteclick to collapse
Some entries in UniProt have one primary (citable) accession number and some secondary numbers, appearing due to merging or demerging entries. In this problem, you may be given any type of ID. If you type the secondary ID into the UniProt query, then you will be automatically redirected to the page containing the primary ID. You can find more information about UniProt IDs here.
"""

from utils.read_fasta import fasta_to_sequence
import re
from typing import Pattern
import requests

# lookahead assertion needed
pattern = re.compile(r"(?=N[^P][ST][^P])")

def match(sequence, pattern: Pattern[str]):
    matches = pattern.finditer(sequence)
    if matches is not None:
        return [m.span()[0] + 1 for m in matches]
    return []

def get_fasta(id: str):
    url = f"https://www.uniprot.org/uniprot/{id}.fasta"
    response = requests.get(url)
    if response.ok:
        return response.text
    else:
        return ""

filename = "./input/rosalind_mprt.txt"
ids = []
with open(filename, "r") as f:
    for line in f.readlines():
        l = line.strip()
        if not l == "":
            ids.append(l)

for id in ids:
    text = get_fasta(id)
    seqs = fasta_to_sequence(text)
    for k, v in seqs.items():
        matches = match(v, pattern)
        if len(matches) > 0:
            print(id)
            print(" ".join(str(m) for m in matches))


