from collections import defaultdict
from typing import DefaultDict, Dict, List, OrderedDict, Tuple
from .recv import reverse_complement



def try_seq_from_dbu(dbu: List[Tuple[str, str]]) -> str:
    build = ""
    # sort it just so outputs are identical
    dbu.sort()
    j = 0
    while len(dbu) > 0 and j < len(dbu):
        curr = dbu[j]
        found = [i for i in dbu if i[0] == curr[1]]
        if found:
            dbu.remove(found[0])
            curr = found[0]
            build += found[0][0][0]
        else:
            j += 1
            continue
    return build

def de_bruin_graph(seqs: List[str]) -> List[Tuple[str, str]]:
    adj_list: List[Tuple[str, str]] = []
    s = set(seqs)
    for seq in seqs:
        s.add(reverse_complement(seq))
    for seq in s:
        adj_list.append((seq[:-1], seq[1:]))
    adj_list.sort()
    return adj_list

def get_kmers(s: str, k: int):
    return [s[i:i +k] for i in range(len(s) - k + 1)]

def find_cycles(reads: List[str], k: int):
    kmers = []
    rc_kmers = []
    for s in reads:
        kmers.extend(get_kmers(s, k))
        kmers.extend(get_kmers(reverse_complement(s), k))
    union_kmers = list(set(kmers) | set(rc_kmers))
    pre_suff = set()
    for s in union_kmers:
        pre_suff.add(s[:-1])
        pre_suff.add(s[1:])
    pre_suff = list(pre_suff)
    nodes: OrderedDict[str, List[str]] = OrderedDict()
    for s in pre_suff:
        nodes[s] = []
    edges: DefaultDict[str, List[str]] = DefaultDict(list)
    for s in union_kmers:
        pre = s[:-1]
        suff = s[1:]
        edges[pre].append(suff)
        edges[suff].append(pre)
    cycles = [i for i in simple_cycles(edges)]
    build = ""
    print(len(cycles))
    print(cycles)
    if len(cycles) == 2:
        print("two")
        print(cycles)
        strings = []
        for c in cycles:
            strings.append(''.join(map(lambda s: s[-1], c)))
        print(strings)

    



def genome_assembly_from_reads(reads: List[str]) -> str:
    for k in range(len(reads[0]), 0, -1):
        find_cycles(reads, k)

    return ""

def simple_cycles(G):
    # Yield every elementary cycle in python graph G exactly once
    # Expects a dictionary mapping from vertices to iterables of vertices
    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    G = {v: set(nbrs) for (v,nbrs) in G.items()} # make a copy of the graph
    sccs = strongly_connected_components(G)
    while sccs:
        scc = sccs.pop()
        startnode = scc.pop()
        path=[startnode]
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [ (startnode,list(G[startnode])) ]
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    stack.append( (nextnode,list(G[nextnode])) )
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            if not nbrs:
                if thisnode in closed:
                    _unblock(thisnode,blocked,B)
                else:
                    for nbr in G[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
        remove_node(G, startnode)
        H = subgraph(G, set(scc))
        sccs.extend(strongly_connected_components(H))

def strongly_connected_components(graph):
    # Tarjan's algorithm for finding SCC's
    # Robert Tarjan. "Depth-first search and linear graph algorithms." SIAM journal on computing. 1972.
    # Code by Dries Verdegem, November 2012
    # Downloaded from http://www.logarithmic.net/pfh/blog/01208083168

    index_counter = [0]
    stack = []
    lowlink = {}
    index = {}
    result = []
    
    def _strong_connect(node):
        index[node] = index_counter[0]
        lowlink[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
    
        successors = graph[node]
        for successor in successors:
            if successor not in index:
                _strong_connect(successor)
                lowlink[node] = min(lowlink[node],lowlink[successor])
            elif successor in stack:
                lowlink[node] = min(lowlink[node],index[successor])

        if lowlink[node] == index[node]:
            connected_component = []

            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            result.append(connected_component[:])
    
    for node in graph:
        if node not in index:
            _strong_connect(node)
    
    return result

def remove_node(G, target):
    # Completely remove a node from the graph
    # Expects values of G to be sets
    del G[target]
    for nbrs in G.values():
        nbrs.discard(target)

def subgraph(G, vertices):
    # Get the subgraph of G induced by set vertices
    # Expects values of G to be sets
    return {v: G[v] & vertices for v in vertices}