from typing import OrderedDict
from fastapi import HTTPException

def fasta_to_sequence(fasta: str) -> OrderedDict[str, str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    last = ""
    for f in fasta.split("\r"):
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
            sequences[last] += line
    return sequences