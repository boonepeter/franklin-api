
from typing import Dict, Tuple, List, Union
from .lcsq import build_longest_matrix
from itertools import product

def edta(s: str, t: str) -> Tuple[int, str, str]:
    sl, tl = len(s), len(t)
    m = {(0, 0): (0, None)}
    m.update({((i, 0), (i, (i-1, 0))) for i in range(1, sl+1)})
    m.update({((0, i), (i, (0, i-1))) for i in range(1, tl+1)})
    i, j = 0, 0
    for i, j in product(range(1, sl+1), range(1, tl+1)):
        cost = 0 if s[i-1] == t[j-1] else 1
        d = m[(i-1, j-1)][0] + cost
        l = m[(i-1, j)][0] + 1
        u = m[(i, j-1)][0] + 1
        b = min(d, l, u)
        if d == b:
            m[(i, j)] = (b, (i-1, j-1))
        elif l == b:
            m[(i, j)] = (b, (i-1, j))
        elif u == b:
            m[(i, j)] = (b, (i, j-1))
    dist = m[(i,j)][0]
    c = (i, j)
    sa = ''
    ta = ''
    while m[c][1] != None:
        i, j = c
        c = m[c][1]
        if i-1 == c[0] and j-1 == c[1]:
            sa += s[i-1]
            ta += t[j-1]
        elif i == c[0]:
            sa += '-'
            ta += t[j-1]
        elif j == c[1]:
            sa += s[i-1]
            ta += '-'
    return dist, ''.join(reversed(sa)), ''.join(reversed(ta))


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