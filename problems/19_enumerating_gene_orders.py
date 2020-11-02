"""
Problem
A permutation of length n is an ordering of the positive integers {1,2,…,n}. For example, π=(5,3,2,1,4) is a permutation of length 5.

Given: A positive integer n≤7.

Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).

Sample Dataset
3
Sample Output
6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1
"""

from itertools import permutations
from typing import List, Tuple


def get_permutations(number: int) -> List[Tuple]:
    perms = permutations(range(1, number + 1))
    return [i for i in perms]


n = 7


with open("./output/permutations.txt", "w") as file:
    perms = get_permutations(n)
    file.write(f"{len(perms)}\n")
    for i in perms:
        file.write(" ".join([str(j) for j in i]) + "\n")

