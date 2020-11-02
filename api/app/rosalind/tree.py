


from typing import List, Tuple

class Node:
    def __init__(self, index: int):
        self.index = index
        self.connections: List[Node] = []
    def __repr__(self) -> str:
        return f"Node({self.index}) - connections: {[i.index for i in self.connections]}"


def build_tree(length: int, nodes: List[Tuple[int, int]]):
    tree = [Node(i) for i in range(1, length + 1)]
    for edge in nodes:
        i = edge[0] - 1
        j = edge[1] - 1
        tree[i].connections.append(tree[j])
        tree[j].connections.append(tree[i])
    return tree
