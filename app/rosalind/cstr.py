

from typing import List

def character_table(t):
    # Return if at the end of the sequence...
    if len(t[0]) < 1:
        return None
        
    # Can't create seperate taxa with only one sequence...
    if len(t) < 2:
        return None
    
    # Initialize the character table.
    char = [[0 for i in t]]
    
    taxa = [[0], []]
    for i in range(1, len(t)):
        if t[i][0] == t[0][0]:
            taxa[0].append(i)
        else:
            taxa[1].append(i)
            
    # Return if all sequences are identical at the current position...
    if taxa[1] == []:
        return None
    
    a, b = [[], []]
    for i in taxa[0]:
        char[0][i] = 1
        a.append(t[i][1:])
    for j in taxa[1]:
        b.append(t[j][1:])
    
    row1 = character_table(a)
    if row1 != None:
        for i in row1:
            char.append(i)
            
    row2 = character_table(b)
    if row2 != None:
        for j in row2:
            char.append(j)
    return char
    


def build_char_table(seqs: List[str]):
    base = seqs[0]
    rows = []
    for i in range(len(base)):
        cnt_0 = 0
        cnt_1 = 0
        row = ''
        for seq in seqs:
            if base[i] == seq[i]:
                row += '1'
                cnt_1 += 1
            else:
                row += '0'
                cnt_0 += 1
        if cnt_1 > 1 and cnt_0 > 1:
            rows.append(row)
    return rows