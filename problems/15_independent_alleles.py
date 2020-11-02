"""
Problem

Figure 2. The probability of each outcome for the sum of the values on two rolled dice (black and white), broken down depending on the number of pips showing on each die. You can verify that 18 of the 36 equally probable possibilities result in an odd sum.
Two events A and B are independent if Pr(A and B) is equal to Pr(A)×Pr(B). In other words, the events do not influence each other, so that we may simply calculate each of the individual probabilities separately and then multiply.

More generally, random variables X and Y are independent if whenever A and B are respective events for X and Y, A and B are independent (i.e., Pr(A and B)=Pr(A)×Pr(B)).

As an example of how helpful independence can be for calculating probabilities, let X and Y represent the numbers showing on two six-sided dice. Intuitively, the number of pips showing on one die should not affect the number showing on the other die. If we want to find the probability that X+Y is odd, then we don't need to draw a tree diagram and consider all possibilities. We simply first note that for X+Y to be odd, either X is even and Y is odd or X is odd and Y is even. In terms of probability, Pr(X+Y is odd)=Pr(X is even and Y is odd)+Pr(X is odd and Y is even). Using independence, this becomes [Pr(X is even)×Pr(Y is odd)]+[Pr(X is odd)×Pr(Y is even)], or (12)2+(12)2=12. You can verify this result in Figure 2, which shows all 36 outcomes for rolling two dice.

Given: Two positive integers k (k≤7) and N (N≤2k). In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.

Return: The probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second law holds for the factors.

Sample Dataset
2 1
Sample Output
0.684
An Example of Dependent Random Variablesclick to collapse

Figure 3. The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is red is equal to (2/5)(1/4), or 1/10.
Two random variables are dependent if they are not independent. For an example of dependent random variables, recall our example in “Mendel's First Law” of drawing two balls from a bag containing 3 red balls and 2 blue balls. If X represents the color of the first ball drawn and Y is the color of the second ball drawn (without replacement), then let A be the event "X is red" and B be the event Y is blue. In this case, the probability tree diagram illustrated in Figure 3 demonstrates that Pr(A and B)=310. Yet we can also see that Pr(A)=35 and Pr(B)=310+110=25. We can now see that Pr(A and B)≠Pr(A)×Pr(B).
"""
from scipy.special import binom

k = 7
N = 32

def foo(k, N):
    def p(n, k):
        return binom(2**k, n) * 0.25**n * 0.75**(2**k - n)
    return 1 - sum(p(n, k) for n in range(N))

print(round(foo(k, N), 3))

