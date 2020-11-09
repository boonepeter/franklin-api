

def probability(length: int, seq: str, x: float) -> float:
    
    gc_count = sum((1 for i in seq if i in "GC"))

    return 1 - (1 - (0.5 * x) ** gc_count * (0.5 * (1 - x)) ** (len(seq) - gc_count)) ** length