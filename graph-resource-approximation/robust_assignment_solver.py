"""Polynomial 4/3 robust graph-resource solver (unit cost, simple ordinary graph).
Arguments: n jobs; edges list of distinct unordered pairs; loop_vertices lists
singleton resources (duplicates allowed); F is a mandatory-resource bit mask.
Resource order is ordinary edges followed by listed singleton resources.
Only singleton bits may be mandatory. Forbid resources by deleting them first.
Returns selected resource bit mask, or None when infeasible.
"""
from polynomial_degree_cover import minimum_degree_cover

def components(n,edges):
    seen=set(); out=[]
    for v in range(n):
        if v in seen: continue
        todo=[v]; seen.add(v); comp=[]
        while todo:
            w=todo.pop(); comp.append(w)
            for a,b in edges:
                u=b if a==w else a if b==w else None
                if u is not None and u not in seen: seen.add(u);todo.append(u)
        out.append(set(comp))
    return out
def solve_graph_resources(n,edges,loop_vertices,F):
    if n<1:raise ValueError('At least one job is required.')
    edges=[tuple(e) for e in edges]
    if any(len(e)!=2 or e[0]==e[1] or any(v<0 or v>=n for v in e) for e in edges):
        raise ValueError('Ordinary edges must have two distinct valid endpoints.')
    if len({frozenset(e) for e in edges})!=len(edges):
        raise ValueError('Parallel ordinary resources are outside this theorem.')
    if any(v<0 or v>=n for v in loop_vertices):raise ValueError('Invalid singleton job.')
    m=len(edges)
    if F<0 or F&((1<<m)-1) or F>>(m+len(loop_vertices)):
        raise ValueError('Mandatory bits must refer only to available singleton resources.')
    res=edges+[(v,) for v in loop_vertices]; answer=0
    for vs in components(n,edges):
        es=[i for i,(a,b) in enumerate(edges) if a in vs]
        ls=[m+j for j,v in enumerate(loop_vertices) if v in vs]
        forced=[i for i in ls if F>>i&1]
        deg={v:sum(v in edges[i] for i in es) for v in vs}
        if len(vs)==1:
            if len(ls)<2:return None
            chosen=set(forced)
            for i in ls:
                if len(chosen)>=2:break
                chosen.add(i)
            answer |=sum(1<<i for i in chosen);continue
        if all(deg[v]==2 for v in vs) and not forced:
            if not ls:return None
            answer |=sum(1<<i for i in es)+(1<<ls[0]);continue
        K=set(forced)
        for v in vs:
            if deg[v]==1 and not any(res[i]==(v,) for i in K):
                opts=[i for i in ls if res[i]==(v,)]
                if not opts:return None
                K.add(opts[0])
        b={v:max(0,2-sum(res[i]==(v,) for i in K)) for v in vs}
        poly_C=minimum_degree_cover(n,[edges[i] for i in es],[b.get(v,0) for v in range(n)])
        if poly_C is None:return None
        C={es[i] for i in poly_C}
        while True:
            bad=[]
            for q in components(n,[edges[i] for i in C]):
                if not q<=vs:continue
                if any(res[i][0] in q for i in K):continue
                if all(sum(v in edges[i] for i in C)==2 for v in q):bad.append(q)
            if not bad:break
            q=bad[0]
            available=[i for i in es if i not in C and any(v in q for v in edges[i])]
            assert available,(n,edges,loop_vertices,F,C,K,q)
            C.add(available[0])
        answer |=sum(1<<i for i in C|K)
    return answer
