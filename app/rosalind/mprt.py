import re
import requests
from typing import Pattern, List
from ..utils.read_fasta import fasta_to_sequence


# lookahead assertion needed
pattern = re.compile(r"(?=N[^P][ST][^P])")

def get_all_motifs(ids: List[str], pattern: Pattern = pattern):
    motifs = {}
    for id in ids:
        text = get_fasta(id)
        seqs = fasta_to_sequence(text)
        for k, v in seqs.items():
            matches = match(v, pattern)
            if len(matches) > 0:
                motifs[id] = " ".join(str(m) for m in matches)

    return motifs



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
