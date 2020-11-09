


from typing import DefaultDict, Tuple, List
from .recv import reverse_complement

def corrections(seqs: List[str]) -> List[Tuple[str, str]]:
    counts = DefaultDict(int)
    for i in seqs:
        counts[i] += 1
        counts[reverse_complement(i)] += 1
    corrs = []
    for i in seqs:
        if counts[i] < 2:
            to_break = False
            for j in range(len(i)):
                if to_break:
                    break
                for bp in "ACGT":
                    potential = i[:j] + bp + i[j + 1:]
                    if potential in counts and counts[potential] > 1:
                        corrs.append((i, potential))
                        to_break = True
                        break
    return corrs