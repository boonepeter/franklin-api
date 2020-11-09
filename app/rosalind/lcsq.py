

from typing import List

def backtrack(C: List[List[int]], seq1: str, seq2: str):
    build = ""
    x = len(seq1)
    y = len(seq2)
    while x > 0 and y > 0:
        if seq1[x - 1] == seq2[y - 1]:
            build = seq1[x - 1] + build
            x -= 1
            y -= 1
            continue
        if C[x][y - 1] > C[x - 1][y]:
            y -= 1
            continue
        x -= 1
    return build

def build_longest_matrix(seq1: str, seq2: str) -> List[List[int]]:
    C = [[0 for i in range(len(seq2) + 1)] for j in range(len(seq1) + 1)]
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            if seq1[i - 1] == seq2[j - 1]:
                C[i][j] = C[i - 1][j-1] + 1
            else:
                C[i][j] = max(C[i-1][j], C[i][j-1])
    return C

def longest_subsequence(seq1: str, seq2: str):
    C = build_longest_matrix(seq1, seq2)
    longest = backtrack(C, seq1, seq2)
    return longest


