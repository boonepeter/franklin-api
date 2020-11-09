


from typing import OrderedDict, Tuple
from ..utils.is_valid import raise_valid_dna

def max_gc_dna(seqs: OrderedDict[str, str]) -> Tuple[str, float]:
    max_gc = 0
    max_name = ""
    for k, v in seqs.items():
        raise_valid_dna(v)
        total = len(v)
        gc = 0
        for l in v:
            if l == "G" or l == "C":
                gc += 1
        percent = gc / total
        if percent > max_gc:
            max_name = k
            max_gc = percent
    return max_name, max_gc