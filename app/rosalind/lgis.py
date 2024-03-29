from typing import List

def subsequence(seq: List[int], increasing: bool=True):
    if not seq:
        return seq

    M = [0] * len(seq)    # offset by 1 (j -> j-1)
    P = [0] * len(seq)

    # Since we have at least one element in our list, we can start by 
    # knowing that the there's at least an increasing subsequence of length one:
    # the first element.
    L = 1
    M[0] = 0

    # Looping over the sequence starting from the second element
    for i in range(1, len(seq)):
        # Binary search: we want the largest j <= L
        #  such that seq[M[j]] < seq[i] (default j = 0),
        #  hence we want the lower bound at the end of the search process.
        lower = 0
        upper = L

        # Since the binary search will not look at the upper bound value,
        # we'll have to check that manually
        if increasing and seq[M[upper-1]] < seq[i]:
                j = upper
        elif not increasing and seq[M[upper - 1]] > seq[i]:
                j = upper
        else:
            # actual binary search loop
            while upper - lower > 1:
                mid = (upper + lower) // 2
                if increasing and seq[M[mid-1]] < seq[i]:
                    lower = mid
                elif not increasing and seq[M[mid-1]] > seq[i]:
                    lower = mid
                else:
                    upper = mid

            j = lower    # this will also set the default value to 0

        P[i] = M[j-1]

        if increasing and (j == L or seq[i] < seq[M[j]]):
            M[j] = i
            L = max(L, j+1)
        elif not increasing and (j == L or seq[i] > seq[M[j]]):
            M[j] = i
            L = max(L, j+1)

    # Building the result: [seq[M[L-1]], seq[P[M[L-1]]], seq[P[P[M[L-1]]]], ...]
    result = []
    pos = M[L-1]
    for _ in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return result[::-1]    # reversing