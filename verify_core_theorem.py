"""Independent, dependency-free verification of bounded robust-core theorem."""
from itertools import combinations
from random import Random
import json

CHECK_DIRECT = False

def pop(x): return x.bit_count()
def matching(rows, right):
    n=len(rows); mate=[-1]*n
    def augment(i,seen):
        for j in range(n):
            if right>>j&1 and rows[i]>>j&1 and not seen[j]:
                seen[j]=True
                if mate[j]<0 or augment(mate[j],seen):
                    mate[j]=i; return True
        return False
    return sum(augment(i,[False]*n) for i in range(n))

def sources(rows):
    n=len(rows); reach=[rows[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    rem=(1<<n)-1; cs=[]
    while rem:
        i=(rem&-rem).bit_length()-1
        c=sum(1<<j for j in range(n) if reach[i]>>j&1 and reach[j]>>i&1)
        rem &= ~c
        if not any(rows[j]&c for j in range(n) if not(c>>j&1)): cs.append(c)
    return cs

def check(rows):
    n=len(rows); v=(1<<n)-1; nu=matching(rows,v); d=n-nu; cs=sources(rows)
    ranks=[pop(s)+matching(rows,v^s)-nu for s in range(v+1)]
    def robust(s): return ranks[s]==d and all(ranks[s^(1<<i)]==d for i in range(n) if s>>i&1)
    def accessible(s): return all(s&c for c in cs)
    def direct_sc(s):
        # Dedicated input rows, each containing one incident state.
        augmented=list(rows)+[1<<i for i in range(n) if s>>i&1]
        # Independent DP over matched right subsets, not dual-rank formula.
        masks={0}
        for row in augmented:
            masks |= {m|1<<j for m in tuple(masks) for j in range(n) if row>>j&1 and not(m>>j&1)}
        return v in masks and accessible(s)
    feasible=[]
    for s in range(v+1):
        robust_full=robust(s) and all(pop(s&c)>=2 for c in cs)
        # Direct independent oracle on sampled graphs, enabled by caller flag.
        if CHECK_DIRECT:
            direct=direct_sc(s) and all(direct_sc(s^(1<<i)) for i in range(n) if s>>i&1)
            assert robust_full==direct,(rows,s,'oracle')
        if robust_full: feasible.append(s)
    opt=min(map(pop,feasible),default=None)
    cores=[s for s in range(v+1) if pop(s)<=2*d and robust(s)]
    predicted=min((pop(k)+sum(max(0,2-pop(k&c)) for c in cs) for k in cores),default=None) if all(pop(c)>=2 for c in cs) else None
    assert opt==predicted,(rows,d,opt,predicted)
    assert (opt is not None)==all(pop(c)>=2 for c in cs),(rows,'feasibility')
    if d==0 and opt is not None: assert opt==2*len(cs)
    if d==1 and opt is not None:
        u=sum(1<<i for i in range(n) if ranks[1<<i])
        q=sum(min(2,pop(u&c)) for c in cs)
        assert opt==2*len(cs)+max(0,2-q),(rows,'d1')
    if d==2:
        minimal=[s for s in cores if not any(robust(s^(1<<i)) for i in range(n) if s>>i&1)]
        for s in minimal:
            ids=[i for i in range(n) if s>>i&1]; classes=[]
            for i in ids:
                assert ranks[1<<i]==1
                for c in classes:
                    if ranks[(1<<i)|(1<<c[0])]==1: c.append(i); break
                else: classes.append([i])
            assert sorted(map(len,classes)) in ([1,1,1],[2,2]),(rows,s,classes)
    return d,opt

if __name__=='__main__':
    counts={}; total=0; CHECK_DIRECT=True
    for n in range(1,4):
        for bits in range(1<<(n*n)):
            rows=[(bits>>(n*i))&((1<<n)-1) for i in range(n)]
            d,opt=check(rows); counts[d]=counts.get(d,0)+1; total+=1
    CHECK_DIRECT=False
    n=4; es=[(i,j) for i in range(n) for j in range(n) if i!=j]
    for bits in range(1<<len(es)):
        rows=[0]*n
        for k,(i,j) in enumerate(es):
            if bits>>k&1: rows[i]|=1<<j
        d,opt=check(rows); counts[d]=counts.get(d,0)+1; total+=1
    rng=Random(20261103)
    for n in range(4,9):
        for trial in range(100):
            p=rng.choice([.12,.25,.5,.75]); rows=[sum(1<<j for j in range(n) if rng.random()<p) for i in range(n)]
            d,opt=check(rows); counts[d]=counts.get(d,0)+1; total+=1
    print(json.dumps({'graphs_checked':total,'deficiency_histogram':counts,'independent_augmented_matching_oracle_graphs':530,'seed':20261103,'result':'PASS'},indent=2))
