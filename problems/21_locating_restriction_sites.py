"""
Problem

Figure 2. Palindromic recognition site
A DNA string is a reverse palindrome if it is equal to its reverse complement. For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC. See Figure 2.

Given: A DNA string of length at most 1 kbp in FASTA format.

Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.

Sample Dataset
>Rosalind_24
TCAATGCATGCGGGTCTATATGCAT
Sample Output
4 6
5 4
6 6
7 4
17 4
18 4
20 6
21 4
Extra Informationclick to collapse
You may be curious how the bacterium prevents its own DNA from being cut by restriction enzymes. The short answer is that it locks itself from being cut through a chemical process called DNA methylation.

"""
from utils.read_fasta import read_fasta

REV = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

def is_palindrom(seq: str) -> bool:
    for i in range(len(seq)):
        a = seq[i]
        b = seq[-(i + 1)]
        if not a == REV[b]:
            return False
    return True

def get_palindrom_locations(sequence: str):
    pals = []
    for i in range(len(sequence)):
        for j in range(4, 13):
            if i + j > len(sequence):
                continue
            s = sequence[i: i + j]
            if is_palindrom(s):
                pals.append((i + 1, len(s)))
    return pals


if __name__ == "__main__":
    seqs = read_fasta("./input/rosalind_revp.txt")
    sequence = [i for i in seqs.values()][0]
    pals = get_palindrom_locations(sequence)
    for i in pals:
        print(" ".join([str(j) for j in i]))