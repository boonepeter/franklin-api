

from typing import Dict


def read_fasta(filename: str) -> Dict[str, str]:
    sequences = {}
    with open(filename, "r") as file:
        last = ""
        for f in file.readlines():
            line = f.strip()
            if f.startswith(">"):
                sequences[line[1:]] = ""
                last = line[1:]
            else:
                sequences[last] += line
    return sequences

def fasta_to_sequence(fasta: str) -> Dict[str, str]:
    sequences = {}
    last = ""
    for f in fasta.split("\n"):
        line = f.strip()
        if line.startswith(">"):
            header = line[1:].split("|")
            name = header[1]
            sequences[name] = ""
            last = name
        else:
            sequences[last] += line
    return sequences