

def fibonacci_rabbits(months: int, pairs: int) -> int:
    """
    Return: The total number of rabbit pairs that will be present after n months,
    if we begin with 1 pair and in each generation, every pair of reproduction-age
    rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
    """
    sequence = [1, 1]
    for i in range(months - 2):
        sequence.append(sequence[i] * pairs + sequence[i + 1])
    return sequence[-1]
