
from typing import List

from fastapi.datastructures import UploadFile


def expected_sites(n: int, seq: str, gc_ratios: List[float]) -> List[float]:
    expected = []

    g = seq.count("G")
    c = seq.count("C")
    t = seq.count("T")
    a = seq.count("A")

    e = n - (len(seq) - 1)
    for i in gc_ratios:
        gc = i * 0.5
        at = (1 - i) * 0.5
        value = (gc ** g) * (gc ** c) * (at ** a) * (at ** t) * e
        expected.append(value)
    return expected

