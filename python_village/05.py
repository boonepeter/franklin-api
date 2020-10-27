"""
Problem
Given: A file containing at most 1000 lines.

Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.

Sample Dataset
Bravely bold Sir Robin rode forth from Camelot
Yes, brave Sir Robin turned about
He was not afraid to die, O brave Sir Robin
And gallantly he chickened out
He was not at all afraid to be killed in nasty ways
Bravely talking to his feet
Brave, brave, brave, brave Sir Robin
He beat a very brave retreat
"""

input_filename = "./inputs/rosalind_ini5.txt"
output_filename = "./output/rosalid_out5.txt"

lines = []

with open(input_filename, "r") as file:
    for line in file.readlines():
        lines.append(line)


with open(output_filename, "w") as file:
    for i in range(len(lines)):
        if (i + 1) % 2 == 0:
            file.write(lines[i])

