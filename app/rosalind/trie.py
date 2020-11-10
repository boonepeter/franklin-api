

from typing import List

class Node:
    Counter = 0
    def __init__(self, value: str) -> None:
        self.value = value
        self.children: List[Node] = []
        Node.Counter += 1
        self.id = Node.Counter
        self.parent: Node = None
    def add_node(self, value: str):
        if len(value) > 0:
            found = False
            for c in self.children:
                if c.value == value[0]:
                    c.add_node(value[1:])
                    found = True
            if not found:
                child = Node(value[0])
                child.parent = self
                self.children.append(child)
                child.add_node(value[1:])
    def __repr__(self) -> str:
        return f"Node({self.value})"
    def print_node(self, line_list):
        if not self.value == "":
            line_list.append(f"{self.parent.id} {self.id} {self.value}")
        for c in self.children:
            line_list.extend(c.print_node([]))
        return line_list

def build_trie(strings: List[str]):
    Node.Counter = 0
    root = Node("")
    for s in strings:
        root.add_node(s)
    lines = []
    root.print_node(lines)
    lines = [i for i in lines]
    return lines
    