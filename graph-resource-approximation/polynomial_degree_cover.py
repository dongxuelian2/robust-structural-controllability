"""Polynomial minimum degree cover via a capacitated matching gadget.
Blossom routine copied unchanged from the independently audited existing
three_halves_solver.py. No third-party packages are required.
"""
from collections import deque

def maximum_matching(adjacency):
    """Return symmetric mate array for a maximum matching of a simple graph.

    Input is a sequence of symmetric neighbor collections without self loops.
    This function does not mutate its input. Each augmenting search contracts
    odd alternating cycles using base representatives, then lifts the path via
    its updated predecessor links. O(v^3) time and O(v^2) input storage.
    """
    n = len(adjacency)
    graph = [tuple(sorted(set(ns))) for ns in adjacency]
    if any(not isinstance(v, int) or v < 0 or v >= n or v == u
           for u, ns in enumerate(graph) for v in ns):
        raise ValueError('Expected simple undirected graph adjacency.')
    neighbor_sets = [set(ns) for ns in graph]
    if any(u not in neighbor_sets[v] for u, ns in enumerate(graph) for v in ns):
        raise ValueError('Expected symmetric adjacency.')
    mate = [-1] * n

    def augment(root):
        base = list(range(n))
        parent = [-1] * n
        outer = [False] * n
        queue = deque([root])
        outer[root] = True

        def ancestor(u, v):
            seen = [False] * n
            while True:
                u = base[u]
                seen[u] = True
                if mate[u] == -1:
                    break
                u = parent[mate[u]]
            while True:
                v = base[v]
                if seen[v]:
                    return v
                v = parent[mate[v]]

        def mark_arm(v, common, child, contracted):
            while base[v] != common:
                contracted[base[v]] = True
                contracted[base[mate[v]]] = True
                parent[v] = child
                child = mate[v]
                v = parent[mate[v]]

        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if base[u] == base[v] or mate[u] == v:
                    continue
                if v == root or (mate[v] != -1 and parent[mate[v]] != -1):
                    common = ancestor(u, v)
                    contracted = [False] * n
                    mark_arm(u, common, v, contracted)
                    mark_arm(v, common, u, contracted)
                    for x in range(n):
                        if contracted[base[x]]:
                            base[x] = common
                            if not outer[x]:
                                outer[x] = True
                                queue.append(x)
                elif parent[v] == -1:
                    parent[v] = u
                    if mate[v] == -1:
                        x = v
                        while x != -1:
                            previous = parent[x]
                            following = mate[previous] if previous != -1 else -1
                            mate[x] = previous
                            if previous != -1:
                                mate[previous] = x
                            x = following
                        return True
                    partner = mate[v]
                    outer[partner] = True
                    queue.append(partner)
        return False

    for root in range(n):
        if mate[root] == -1:
            augment(root)
    return tuple(mate)



def minimum_degree_cover(n,edges,requirements):
    """Return ordinary edge indices, or None if infeasible; O((n+m)^3)."""
    degree=[0]*n
    for u,v in edges:degree[u]+=1;degree[v]+=1
    cap=[degree[v]-requirements[v] for v in range(n)]
    if min(cap,default=0)<0:return None
    copies=[];count=0
    for c in cap:
        copies.append(list(range(count,count+c)));count+=c
    ports=[(count+2*i,count+2*i+1) for i in range(len(edges))]
    graph=[set() for _ in range(count+2*len(edges))]
    def link(a,b):graph[a].add(b);graph[b].add(a)
    for (u,v),(a,b) in zip(edges,ports):
        link(a,b)
        for w in copies[u]:link(a,w)
        for w in copies[v]:link(b,w)
    mate=maximum_matching(graph)
    excluded={i for i,((u,v),(a,b)) in enumerate(zip(edges,ports))
              if mate[a] in copies[u] and mate[b] in copies[v]}
    assert len(excluded)==sum(x>=0 for x in mate)//2-len(edges)
    result=set(range(len(edges)))-excluded
    assert all(sum(v in edges[i] for i in result)>=requirements[v] for v in range(n))
    return result
