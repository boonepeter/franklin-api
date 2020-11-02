

from itertools import permutations, product

filename = "./input/rosalind_sign.txt"

num = 5


with open(filename, "r") as f:
    contents = f.read().strip()
    num = int(contents)

pos = [i for i in range(1, num + 1)]
neg = [-i for i in range(1, num + 1)]

prod = product(pos, neg, repeat=num)



nums = [(i, -i) for i in range(1, num + 1)]
nums = [i for sublist in nums for i in sublist]
perms = list(permutations(nums, r=num))

good = []
for p in perms:
    bad = False
    for i in range(1, num + 1):
        if i in p and -i in p:
            bad = True
            break
    if not bad:
        good.append(p)

for g in good:
    print(g)
resp = f"{len(good)}\n"
resp += "\n".join(" ".join(str(j) for j in i) for i in good)

with open("./output/29.txt", "w") as out:
    out.write(resp)
