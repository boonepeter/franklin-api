


from typing import List, Tuple

from .recv import reverse_complement

def de_bruin_graph(seqs: List[str]) -> List[Tuple[str, str]]:
    adj_list: List[Tuple[str, str]] = []
    s = set(seqs)
    for seq in seqs:
        s.add(reverse_complement(seq))
    for seq in s:
        adj_list.append((seq[:-1], seq[1:]))
    adj_list.sort()
    return adj_list