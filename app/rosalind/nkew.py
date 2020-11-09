
from .nwck import Node, Parser
from typing import List



def pathToNode(root: Node, path: List[Node], k: str): 
    # base case handling
    if root is None:
        return False
 
     # append the node value in path
    path.append(root)
  
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

def path_length(nodes: List[Node]) -> int:
    total = 0
    for i in nodes:
        total += i.branch_length
    return total

def distance(root: Node, data1: str, data2: str) -> int:
    if root:
        # store path corresponding to node: data1
        path1: List[Node] = []
        pathToNode(root, path1, data1)
 
        # store path corresponding to node: data2
        path2: List[Node] = []
        pathToNode(root, path2, data2)
 
        # iterate through the paths to find the 
        # common path length
        i = 0
        dist = 0
        while i < len(path1) and i < len(path2):
            # get out as soon as the path differs 
            # or any path's length get exhausted
            if path1[i] != path2[i]:
                break
            dist += path1[i].branch_length
            i += 1
 
        # get the path length by deducting the 
        # intersecting path length (or till LCA)
        return (path_length(path1) + path_length(path2) - 2 * dist)
    else:
        return 0


def parse_return_distance(input_str: str, node1: str, node2: str) -> int:
    parser = Parser(input_str)
    root =  parser.parse()
    dist = distance(root, node1, node2)
    return dist

