"""
Problem

Figure 4. A figure illustrating the propagation of Fibonacci's rabbits if they die after three months.
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2 and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n≤100 and m≤20.

Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.

Sample Dataset
6 3
Sample Output
4
"""

def mortal_fib(n: int, m: int) -> int:
    # Create an array of length m to keep track of the number of rabbits
    ages = [0] * m
    # Start with one rabbit being alive with m months to live
    ages[-1] = 1
    # Iterate over the number of n months to track
    for _ in range(1, n):
        # Newborns
        new_rabbits = sum(ages[:-1])
        # Shift ages left, i.e. getting older
        ages[:-1] = ages[1:] 
        # Assign newborns 
        ages[-1] = new_rabbits
    
    # At the end we have an array containing the number of rabbits alive
    # for each possible month in the rabbits lifespan (m).
    # The sum of all ages is then a representation of all rabbits
    # currently alive.
    return sum(ages)

if __name__ == '__main__':
    print(mortal_fib(80, 16))  # Prints 4
