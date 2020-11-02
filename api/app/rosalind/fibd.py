

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
