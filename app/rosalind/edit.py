

def levenshtein_distance(seq1: str, seq2: str) -> int:
    """Wagner–Fischer algorithm
    """
    dists = [[0 for i in range(len(seq2))] for j in range(len(seq1))]
    for i in range(len(seq1)):
        dists[i][0] = i
    
    for j in range(len(seq2)):
        dists[0][j] = j

    for j in range(1, len(seq2)):
        for i in range(1, len(seq1)):
            if seq1[i - 1] == seq2[j - 1]:
                sub_cost = 0
            else:
                sub_cost = 1
            
            deletion = dists[i-1][j] + 1
            insertion = dists[i][j -1] + 1
            substitution = dists[i-1][j-1] + sub_cost
            dists[i][j] = min(deletion, substitution, insertion)
    return dists[len(seq1) - 1][len(seq2) - 1]
