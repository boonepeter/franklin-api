from app.utils.is_valid import raise_valid_dna
from typing import OrderedDict
from fastapi import HTTPException

def fasta_to_sequence(fasta: str, raise_dna: bool=True) -> OrderedDict[str, str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    last = ""
    for f in fasta.split("\n"):
        line = f.strip()
        if line.startswith(">"):
            header = line[1:].split("|")
            if len(header) > 1:
                name = header[1]
            else:
                name = line[1:]
            sequences[name] = ""
            last = name
        else:
            if last in sequences:
                sequences[last] += line
            else:
                raise HTTPException(status_code=500, detail="Invalid Fasta format")
    return sequences