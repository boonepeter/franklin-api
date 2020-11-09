


from typing import List


def kmp_failure(seq: str) -> List[int]:
    ret = [0]
    for i in range(1, len(seq)):
        j = ret[i - 1]
        while j > 0 and seq[j] != seq[i]:
            j = ret[j - 1]
        ret.append(j + 1 if seq[j] == seq[i] else j)
    return ret

