
from typing import Dict

def helper(seq: str, memo: Dict[str, int], mod: int=10**6):
    if len(seq) == 0 or len(seq) == 1:
        return 1
    
    total = 0
    if seq in memo:
        return memo[seq]
    
    for i in range(1, len(seq), 2):
        if ((seq[0] == 'A' and seq[i] == 'U') or
            (seq[0] == 'U' and seq[i] == 'A') or
            (seq[0] == 'C' and seq[i] == 'G') or
            (seq[0] == 'G' and seq[i] == 'C')):
            total += helper(seq[1:i], memo, mod) * helper(seq[i+1:], memo, mod)

    memo[seq] = total % 10**6
    return memo[seq]
    
def catalan(seq: str, mod: int=10**6):
    memo = {}
    return helper(seq, memo, mod)
