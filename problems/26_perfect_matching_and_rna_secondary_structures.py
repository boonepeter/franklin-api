"""
Problem

Figure 2. Three matchings (highlighted in red) shown in three different graphs.

Figure 3. A graph containing 10 nodes; the five edges forming a perfect matching on these nodes are highlighted in red.

Figure 4. The bonding graph for the RNA string s = UAGCGUGAUCAC.

Figure 5. A perfect matching on the basepair edges is highlighted in red and represents a candidate secondary structure for the RNA strand.
A matching in a graph G is a collection of edges of G for which no node belongs to more than one edge in the collection. See Figure 2 for examples of matchings. If G contains an even number of nodes (say 2n), then a matching on G is perfect if it contains n edges, which is clearly the maximum possible. An example of a graph containing a perfect matching is shown in Figure 3.

First, let Kn denote the complete graph on 2n labeled nodes, in which every node is connected to every other node with an edge, and let pn denote the total number of perfect matchings in Kn. For a given node x, there are 2n−1 ways to join x to the other nodes in the graph, after which point we must form a perfect matching on the remaining 2n−2 nodes. This reasoning provides us with the recurrence relation pn=(2n−1)⋅pn−1; using the fact that p1 is 1, this recurrence relation implies the closed equation pn=(2n−1)(2n−3)(2n−5)⋯(3)(1).

Given an RNA string s=s1…sn, a bonding graph for s is formed as follows. First, assign each symbol of s to a node, and arrange these nodes in order around a circle, connecting them with edges called adjacency edges. Second, form all possible edges {A, U} and {C, G}, called basepair edges; we will represent basepair edges with dashed edges, as illustrated by the bonding graph in Figure 4.

Note that a matching contained in the basepair edges will represent one possibility for base pairing interactions in s, as shown in Figure 5. For such a matching to exist, s must have the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

Given: An RNA string s of length at most 80 bp having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

Return: The total possible number of perfect matchings of basepair edges in the bonding graph of s.

Sample Dataset
>Rosalind_23
AGCUAGUCAU
Sample Output
12
"""
from math import factorial
rna = "GCAUCGAAUCUCAUUACAUAUUGAGCCGCGCCUCUGUUGCUCAAGGGCAGGAGCCCGCAAUUUUCAGCAGGAGG"
pmch = factorial(rna.count('A'))*factorial(rna.count('C'))

print(pmch)