"""
Problem
For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.

By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).

The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).

Sample Dataset
>Rosalind_56
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC
Sample Output
ATTAGACCTGCCGGAATAC
"""

from typing import List
from utils.read_fasta import read_fasta

fasta = ["ATTAGACCTG",
"CCTGCCGGAA",
"AGACCTGCCG",
"GCCGGAATAC"]

fn = "./input/rosalind_long.txt"

seqs = read_fasta(fn).values()

seqs = list(seqs)


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

if __name__ == "__main__":
    print(find_overlaps(seqs))