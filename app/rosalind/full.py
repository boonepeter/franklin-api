

import re
from typing import List
from .prtm import masses


m_to_p = { round(v, 5): k for k, v in masses.items()}


def infer_from_full(parent: float, ions: List[float]):
    protein = ""
    length = (len(ions) - 2) // 2
    new_i = 0
    while len(protein) < length:
        to_break = False
        for i in range(new_i, len(ions) - 1):
            if to_break:
                break
            for j in range(i, len(ions)):
                dif = round(ions[j] - ions[i], 5)
                if dif in m_to_p:
                    protein += m_to_p[dif]
                    new_i = j
                    to_break = True
                    break
    return protein


