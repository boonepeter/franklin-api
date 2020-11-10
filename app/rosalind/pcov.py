

from typing import List, Tuple


def perfect_coverage(seqs: List[str]) -> str:
    seq_set = set(seqs)
    adj_list: List[Tuple[str, str]] = []
    for i in seq_set:
        adj_list.append((i[:-1], i[1:]))
    build = ""
    count = 0
    # sort it just so outputs are identical
    adj_list.sort()
    curr = adj_list[0]
    while len(adj_list) > 0:
        count += 1
        found = [i for i in adj_list if i[0] == curr[1]]
        if found:
            adj_list.remove(found[0])
            curr = found[0]
            build += found[0][0][0]
        else:
            break
    return build