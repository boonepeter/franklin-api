

from typing import Dict

from collections import OrderedDict


def read_fasta(filename: str):
    sequences = OrderedDict()
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

def fasta_to_sequence(fasta: str):
    sequences = OrderedDict()
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