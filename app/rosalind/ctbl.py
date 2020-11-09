
import enum
from typing import List
from .nwck import get_all_nodes, parse_newick


NODES = []

class Node:
    def __init__(self, par):
        self.p = par
        self.s = set()
        self.lab = ''
        self.sonlabs = set()
        NODES.append(self)

    def __repr__(self):
        return '%s(s=%s)' % (self.lab or 'Node', self.s)


def buildtree(tree):
    class cur:
        pos = 0
    dd = {}

    def getnode(par):
        cc = Node(par)
        if tree[cur.pos] == '(':
            while tree[cur.pos] in '(,':
                cur.pos += 1
                cc.s.add(getnode(cc))
            cur.pos += 1
        ff = cur.pos
        while tree[cur.pos] not in '), ;':
            cur.pos += 1
        nam = tree[ff:cur.pos]
        cc.lab = nam
        if nam != '':
            dd[nam] = cc
        return cc
    return (getnode(None), dd)


def cnt(cur):
    trt = set()
    for son in cur.s:
        trt.update(cnt(son))
    if cur.lab:
        trt.add(cur.lab)
    cur.sonlabs = trt
    return trt


def build_character_table(in_str: str) -> List[List[int]]:
    root = parse_newick(in_str)
    nodes = get_all_nodes(root)
    nodes = [n.name for n in nodes]
    nodes.sort()

    (root, d) = buildtree(in_str)
    alllabs = sorted(cnt(root))
    nn = len(alllabs)
    print(NODES)
    for j in NODES:
        if len(j.sonlabs) not in (0, 1, nn - 1, nn):
            print(''.join(map(str, map(int, map(j.sonlabs.__contains__, alllabs)))))

    return [[]]
