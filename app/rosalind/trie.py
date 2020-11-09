

from typing import List

class Tree:
    def __init__(self) -> None:
        self.root = None

def build_trie(strings: List[str]):
    n = max(len(i) for i in strings)
    lists = []