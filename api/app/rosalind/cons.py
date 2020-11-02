
from collections import defaultdict
from typing import OrderedDict

def find_consensus(fasta: OrderedDict[str, str]) -> str:
    seqs = list(fasta.values())
    consensus = ""
    cons_matrix = defaultdict(list)
    for i in range(len(seqs[0])):
        counts = {
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0
        }
        for j in seqs:
            counts[j[i]] += 1
        max_l = ""
        max_v = -1
        for k, v in counts.items():
            if v > max_v:
                max_l = k
                max_v = v
            cons_matrix[k].append(v)
        consensus += max_l
    #consensus += "\n"
    #for k in "ACGT":
    #    consensus += k + ": " + " ".join(str(i) for i in cons_matrix[k]) + "\n"

    return consensus