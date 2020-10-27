"""
Problem
Given: Two positive integers a and b (a<b<10000).

Return: The sum of all odd integers from a through b, inclusively.
"""

a = 4729
b = 8738

total = 0
for i in range(a, b + 1):
    if i % 2 != 0:
        total += i

print(total)
