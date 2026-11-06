"""Independent finite verification. Standard Python; exact enumeration replaces b-matching."""
import itertools, random
from polynomial_degree_cover import minimum_degree_cover
from robust_assignment_solver import solve_graph_resources

def covers(n,res,mask):
    reach={0}
    for i,ends in enumerate(res):
        if mask>>i&1:
            reach |= {s|1<<v for s in tuple(reach) for v in ends if not s>>v&1}
    return (1<<n)-1 in reach

def robust(n,res,mask):
    return bool(mask) and all(covers(n,res,mask^(1<<i)) for i in range(len(res)) if mask>>i&1)

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

def characterization(n,res,mask):
    deg=[0]*n; edges=[]; loops=[]
    for i,ends in enumerate(res):
        if mask>>i&1:
            for v in ends:deg[v]+=1
            if len(ends)==2:edges.append(ends)
            else:loops.append(ends[0])
    if min(deg)<2:return False
    for vs in components(n,edges):
        if not any(v in vs for v in loops) and all(deg[v]==2 for v in vs):return False
    return True

def algorithm(n,edges,loop_vertices,F):
    m=len(edges);res=edges+[(v,) for v in loop_vertices]; answer=0
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
        candidates=[]
        for x in range(1<<len(es)):
            sel={es[j] for j in range(len(es)) if x>>j&1}
            if all(sum(v in edges[i] for i in sel)>=b[v] for v in vs):candidates.append(sel)
        if not candidates:return None
        exact_C=min(candidates,key=lambda x:(len(x),sorted(x)))
        poly_C=minimum_degree_cover(n,[edges[i] for i in es],[b.get(v,0) for v in range(n)])
        assert poly_C is not None and len(poly_C)==len(exact_C)
        C={es[i] for i in poly_C}
        # Independently verify the complement capacitated-b-matching identity.
        max_excluded=max(sum(x>>j&1 for j in range(len(es)))
            for x in range(1<<len(es))
            if all(sum(bool(x>>j&1) and v in edges[es[j]] for j in range(len(es)))
                   <=deg[v]-b[v] for v in vs))
        assert len(C)==len(es)-max_excluded
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

def fixture(n,edges,loops,forced_patterns):
    res=edges+[(v,) for v in loops]; good=[]
    for x in range(1<<len(res)):
        r=robust(n,res,x)
        assert r==characterization(n,res,x),(n,res,x,r)
        if r:good.append(x)
    tested=0
    for F in forced_patterns:
        feasible=[x for x in good if x&F==F]
        A=algorithm(n,edges,loops,F)
        assert A==solve_graph_resources(n,edges,loops,F)
        if not feasible:assert A is None
        else:
            opt=min(x.bit_count() for x in feasible)
            assert A is not None and A&F==F and robust(n,res,A)
            loopmask=((1<<len(res))-1)^((1<<len(edges))-1)
            canonical=A&loopmask
            assert min(x.bit_count() for x in feasible if x&loopmask==canonical)==opt
            assert 3*A.bit_count()<=4*opt,(n,edges,loops,F,A,opt)
        tested+=1
    return tested,1<<len(res)

def main():
    count=0; masks=0
    for n in range(1,5):
        full=list(itertools.combinations(range(n),2))
        for g in range(1<<len(full)):
            edges=[e for j,e in enumerate(full) if g>>j&1]
            loops=list(range(n)); m=len(edges)
            c,k=fixture(n,edges,loops,[f<<m for f in range(1<<n)])
            count+=c;masks+=k
    rng=random.Random(438519)
    for _ in range(180):
        n=rng.randint(1,5);edges=[e for e in itertools.combinations(range(n),2) if rng.random()<.45]
        loops=[v for v in range(n) for _ in range(rng.randint(0,2))]
        if len(edges)+len(loops)>14:continue
        m=len(edges)
        fs={0}|{sum(1<<(m+j) for j in range(len(loops)) if rng.random()<.35) for _ in range(8)}
        c,k=fixture(n,edges,loops,fs);count+=c;masks+=k
    print(f'PASS: {count} mandatory-loop fixtures; {masks} selected resource subsets independently matched and compared with characterization.')
if __name__=='__main__':main()
