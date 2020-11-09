

from .lcsq import longest_subsequence

def shortest_common_supersequence(seq1: str, seq2: str) -> str:
    build = ""
    lcs = longest_subsequence(seq1, seq2)
    i, j = 0, 0
    for l in lcs:
        while not seq1[i] == l:
            build += seq1[i]
            i += 1
        while not seq2[j] == l:
            build += seq2[j]
            j += 1
        build += l
        j += 1
        i += 1
    build += seq1[i:]
    build += seq2[j:]
    return build