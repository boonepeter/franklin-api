
from typing import OrderedDict


def longest_common_substring(fastas: OrderedDict[str, str]) -> str:
    seqs = list(fastas.values())
    subs = generate_subs(seqs[0])
    longest = 0
    longest_substring = ""

    for s in subs:
        in_all = True
        for f in seqs:
            if s not in f:
                in_all = False
                break
        if in_all:
            if len(s) > longest:
                longest = len(s)
                longest_substring = s
    return longest_substring


def generate_subs(word):
    subs = set()
    for i in range(len(word)):
        for j in range(i + 1, len(word)):
            s = word[i: j]
            subs.add(s)
    return subs