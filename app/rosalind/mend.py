

from .nwck import Parser, Node
from typing import List


def create_probs(current: Node):
    for n in current.children:
        create_probs(n)
    if current.name:
        for i in ["AA", "Aa", "aa"]:
            if i == ''.join(sorted(current.name)):
                current.probs[i] = 1.0
            else:
                current.probs[i] = 0.0
    else:
        # we are assuming that there are no unlabeled leaves.
        # Otherwise we may not be able to infer probability
        c1 = current.children[0].probs
        c2 = current.children[1].probs
        AA = c1['AA'] * c2['AA']
        AA += c1["Aa"] * c2["Aa"] / 4
        AA += c1['AA'] * c2["Aa"] / 2
        AA += c1['Aa'] * c2["AA"] / 2
        Aa = c1['AA'] * c2['aa'] 
        Aa += c1['aa'] * c2['AA'] 
        Aa += c1['Aa'] * c2['Aa'] / 2
        Aa += c1['AA'] * c2['Aa'] / 2
        Aa += c1['Aa'] * c2['AA'] / 2 
        Aa += c1['aa'] * c2['Aa'] / 2 
        Aa += c1['Aa'] * c2['aa'] / 2
        aa = c1['aa'] * c2['aa']
        aa += c1['Aa'] * c2['Aa'] / 4
        aa += c1['aa'] * c2['Aa'] / 2
        aa += c1['Aa'] * c2['aa'] / 2
        current.probs["AA"] = AA
        current.probs["Aa"] = Aa
        current.probs["aa"] = aa


def parse_get_probs(in_str: str) -> List[float]:
    parser = Parser(in_str, auto_name=False)
    root = parser.parse()
    create_probs(root)
    return [root.probs["AA"], root.probs["Aa"], root.probs["aa"]]