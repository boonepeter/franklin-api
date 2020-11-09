

from typing import List

def find_overlaps(arr: List[str], acc: str='') -> str:
    # if no sequences left, return the accumulated string
    if len(arr) == 0:
        return acc

    # accumulated is "", add the first sequence and then go down recursively
    elif len(acc) == 0:
        acc = arr.pop(0)
        return find_overlaps(arr, acc)

    else:
        for i in range(len(arr)):
            a = arr[i]
            l = len(a)

            # we know they overlap by at least half
            for p in range(l // 2):
                q = l - p

                # if accumulated string starts with the last half of a,
                # add a to beginning of acc
                if acc.startswith(a[p:]):
                    arr.pop(i)
                    return find_overlaps(arr, a[:p] + acc)

                if acc.endswith(a[:q]):
                    arr.pop(i)
                    return find_overlaps(arr, acc + a[q:])
    return ""