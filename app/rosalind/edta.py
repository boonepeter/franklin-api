
from typing import Dict, Tuple, List, Union
from .lcsq import build_longest_matrix



def build_inserts(C: List[List[int]], seq1: str, seq2: str):
    build = ""
    build1 = ""
    build2 = ""
    x = len(seq1)
    y = len(seq2)
    count = 0
    while x > 0 and y > 0:
        if seq1[x - 1] == seq2[y - 1]:
            build1 = seq1[x - 1] + build1
            build2 = seq1[x - 1] + build2
            x -= 1
            y -= 1
            continue
        if C[x][y - 1] > C[x - 1][y]:
            build1 = "-" + build1
            build2 = seq2[y - 1] + build2
            y -= 1
        else:
            build2 = "-" + build2
            build1 = seq1[x - 1] + build1
            x -= 1
        count += 1
    return build1, build2, count


def min_edit_distance(seq1: str, seq2: str) -> Tuple[int, str, str]:
    matrix = build_longest_matrix(seq1, seq2)
    b1, b2, c = build_inserts(matrix, seq1, seq2)
    return c, b1, b2