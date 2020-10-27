"""
Problem
A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

Sample Dataset
>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA
Sample Output
AC
"""
# this is a bad solution. Very slow. Need to use suffix trees and stuff. more work needed.

from utils.read_fasta import read_fasta

filename = "./input/rosalind_lcsm.txt"


def longest_common_substring(filename):
    fastas = read_fasta(filename)

    fastas = list(fastas.values())
    subs = generate_subs(fastas[0])
    longest = 0
    longest_substring = ""

    for s in subs:
        in_all = True
        for f in fastas:
            if s not in f:
                in_all = False
                break
        if in_all:
            if len(s) > longest:
                longest = len(s)
                longest_substring = s
    return longest_substring


def generate_subs(word):
    subs = set()
    for i in range(len(word)):
        for j in range(i + 1, len(word)):
            s = word[i: j]
            subs.add(s)
    return subs

print(longest_common_substring(filename))