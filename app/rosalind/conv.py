

from typing import Set, Tuple


def convolution(set1: Set[float], set2: Set[float]) -> Tuple[int, float]:
    ll = sorted(i - j for i in set1 for j in set2)
    mxc = 0
    mxr = 0
    t = 0
    prev = -10000000000000.0
    for i in ll:
        if i - prev > 0.0000001:
            t = 0
        t += 1
        prev = i
        if mxc < t:
            mxc = t
            mxr = i
    return mxc, mxr