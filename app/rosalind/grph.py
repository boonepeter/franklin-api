



from typing import List, OrderedDict


def find_adjacent(fastas: OrderedDict[str, str], dist=3) -> List[str]:
    graphs = []
    for k, v in fastas.items():
        for k2, v2 in fastas.items():
            if k == k2:
                continue
            if v[-dist:] == v2[:dist]:
                graphs.append(k + " " + k2)
    return graphs