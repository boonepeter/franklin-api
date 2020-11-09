import re
from typing import Dict, List, Union

tokens = [
    (r"\(", "open parens"),
    (r"\)", "close parens"),
    (r"[^\s\(\)\[\]\'\:\;\,]+", "unquoted node label"),
    (r"\:\ ?[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?", "edge length"),
    (r"\,", "comma"),
    (r"\[(\\.|[^\]])*\]", "comment"),
    (r"\'(\\.|[^\'])*\'", "quoted node label"),
    (r"\;", "semicolon"),
    (r"\n", "newline"),
]

tokenizer = re.compile("(%s)" % "|".join(t[0] for t in tokens))

token_dict = {name: re.compile(t) for t, name in tokens}


class Node:
    def __init__(self, name: str) -> None:
        self.name = name
        self.children: List[Node] = []
        self.parent: Union[Node, None] = None
        self.branch_length = 0
        self.probs: Dict[str, float] = {}
    def __repr__(self) -> str:
        return f"Node({self.name}, {len(self.children)} children)"
    def print_node(self):
        print(self)
        for c in self.children:
            c.print_node()

def parse_newick(input_string) -> Node:
    parser = Parser(input_string)
    root = parser.parse()
    return root

def get_all_nodes(root: Node) -> List[Node]:
    nodes = []
    nodes.append(root)
    for c in root.children:
        nodes.extend(get_all_nodes(c))
    return nodes



def pathToNode(root: Node, path: List[str], k: str): 
    # base case handling
    if root is None:
        return False
 
     # append the node value in path
    path.append(root.name)
  
    # See if the k is same as root's data
    if root.name == k:
        return True
  
    # Check if k is found in left or right 
    # sub-tree
    if any(pathToNode(c, path, k) for c in root.children):
        return True
  
    # If not present in subtree rooted with root, 
    # remove root from path and return False 
    path.pop()
    return False

def distance(root: Node, data1: str, data2: str):
    if root:
        # store path corresponding to node: data1
        path1 = []
        pathToNode(root, path1, data1)
 
        # store path corresponding to node: data2
        path2 = []
        pathToNode(root, path2, data2)
 
        # iterate through the paths to find the 
        # common path length
        i=0
        while i < len(path1) and i < len(path2):
            # get out as soon as the path differs 
            # or any path's length get exhausted
            if path1[i] != path2[i]:
                break
            i = i + 1
 
        # get the path length by deducting the 
        # intersecting path length (or till LCA)
        return (len(path1) + len(path2) - 2 * i)
    else:
        return 0


class Parser:
    """
    Parse a newick tree format given a string
    """

    def __init__(self, input_string: str, auto_name: bool=True):
        self.contents = input_string
        self.counter = 0
        self.auto_name = auto_name

    def parse(self):
        tokens = re.finditer(tokenizer, self.contents)
        root = self.new_node()
        current = root
        entering = False
        lp_count = 0
        rp_count = 0
        for match in tokens:
            token = match.group()
            if token.startswith("'"):
                current.name = token[1:-1]
            
            elif token.startswith("["):
                # comment
                pass
            elif token == "(":
                current = self.new_node(current)
                entering = False
                lp_count += 1
            
            elif token == ",":
                if current is root:
                    root = self.new_node()
                    current.parent = root
                parent = self.process_node(current)
                current = self.new_node(parent)
                entering = False
            elif token == ")":
                parent = self.process_node(current)
                if not parent:
                    raise Exception("Par missing")
                current = parent
                entering = False
                rp_count += 1
            elif token == ";":
                break
            elif token.startswith(":"):
                value = int(token[1:])
                current.branch_length = value
            elif token == "\n":
                pass
            else:
                current.name = token
        
        if not lp_count == rp_count:
            raise Exception("Par counts don't match")
        self.process_node(current)
        self.process_node(root)
        return root
    def new_node(self, parent: Union[Node, None] = None) -> Node:
        name = str(self.counter) if self.auto_name else "";
        node = Node(name)
        if parent:
            node.parent = parent
        self.counter += 1
        return node
    
    def process_node(self, node: Node):
        try:
            parent = node.parent
        except AttributeError:
            pass
        else:
            if parent:
                parent.children.append(node)
            del node.parent
            return parent

        