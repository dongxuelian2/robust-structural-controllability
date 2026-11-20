"""Four-source-state control-graph checks. No external dependencies."""
from itertools import combinations
from random import Random
from verify_four_thirds import algorithm, robust, components


def matching_size(rows,n):
    mate=[-1]*n
    def aug(i,seen):
        for v in range(n):
            if rows[i]>>v&1 and v not in seen:
                seen.add(v)
                if mate[v]<0 or aug(mate[v],seen):mate[v]=i;return True
        return False
    return sum(aug(i,set()) for i in range(len(rows)))


def build(n,edges):
    resources=[(1<<u)|(1<<v) for u,v in edges]+[1,1<<(n-1)]
    m=len(resources);A,B,FA,FB=range(2*m+n,2*m+n+4);rows=[0]*(2*m+n+4)
    for e,nb in enumerate(resources):
        rows[e]=1<<(m+e);rows[m+e]=(1<<e)|(nb<<(2*m))
    rows[A]=(1<<FA)|(((1<<m)-1)<<m)
    rows[FA]=(1<<B)|(1<<(m-1));rows[B]=1<<FB
    rows[FB]=(1<<A)|(1<<(m-2))
    return rows,(A,B,FA,FB),m


def actual(rows,selected):
    n=len(rows)
    for failed in [-1]+list(selected):
        survivors=set(selected)-{failed}; reach=set(survivors);todo=list(reach)
        while todo:
            v=todo.pop()
            for u in range(n):
                if rows[v]>>u&1 and u not in reach:reach.add(u);todo.append(u)
        if len(reach)!=n:return False
        if matching_size(rows+[1<<v for v in survivors],n)!=n:return False
    return True


def fixture(n,edges,rng):
    rows,source,m=build(n,edges);A,B,FA,FB=source
    # Resource IDs and physical state IDs are deliberately kept separate.
    loops=[0,n-1]+list(range(n))+[0,n-1]
    loop_states=[len(edges),len(edges)+1]+list(range(2*m,2*m+n))+[A,B]
    states=list(range(len(edges)))+loop_states
    res=edges+[(v,) for v in loops]
    ai=len(res)-2;bi=len(res)-1
    good=[x for x in range(1<<len(res)) if robust(n,res,x)]
    best_opt=10**9;best_sol=None;patterns=0
    for z in range(16):
        if z.bit_count()<2:continue
        mandatory=((z&1)>0)<<ai | ((z&2)>0)<<bi
        forbidden=((z&1)==0)<<ai | ((z&2)==0)<<bi
        compatible=[x for x in good if x&mandatory==mandatory and not x&forbidden]
        if not compatible:continue
        charge=bool(z&4)+bool(z&8)
        exact=min(x.bit_count() for x in compatible)+charge
        best_opt=min(best_opt,exact)
        available=[j for j in range(len(loops)) if not (j==len(loops)-2 and not z&1) and not (j==len(loops)-1 and not z&2)]
        currloops=[loops[j] for j in available]
        F=sum(1<<(len(edges)+k) for k,j in enumerate(available) if j>=len(loops)-2)
        ans=algorithm(n,edges,currloops,F)
        assert ans is not None
        currstates=list(range(len(edges)))+[loop_states[j] for j in available]
        selected={currstates[i] for i in range(len(currstates)) if ans>>i&1}
        if z&4:selected.add(FA)
        if z&8:selected.add(FB)
        assert set(source)&selected=={source[i] for i in range(4) if z>>i&1}
        assert actual(rows,selected),(n,edges,z,selected)
        assert 3*len(selected)<=4*exact
        if best_sol is None or len(selected)<len(best_sol):best_sol=selected
        patterns+=1
    assert 3*len(best_sol)<=4*best_opt
    # Random placements include centers, all four source states, and arbitrary sizes.
    checks=0
    for _ in range(140):
        selected={v for v in range(len(rows)) if rng.random()<.4}
        x=sum(1<<i for i,v in enumerate(states) if v in selected)
        expected=len(selected&set(source))>=2 and robust(n,res,x)
        assert actual(rows,selected)==expected
        checks+=1
    # Exhaust actual state placements up to the exact optimum on tiny instances.
    if n<=3:
        found=None
        for k in range(n+1,best_opt+1):
            for selected in combinations(range(len(rows)),k):
                if len(set(selected)&set(source))<2:continue
                value=actual(rows,selected);checks+=1
                if value:found=k;break
            if found is not None:break
        assert found==best_opt
    return patterns,checks,len(good)


def main():
    rng=Random(765181);fixtures=[]
    for n in range(2,5):
        full=list(combinations(range(n),2))
        for mask in range(1<<len(full)):
            edges=[e for i,e in enumerate(full) if mask>>i&1]
            deg=[sum(v in e for e in edges) for v in range(n)]
            if min(deg)==0:continue
            if n==4 and min(deg)<2:continue
            fixtures.append((n,edges))
    fixtures.extend([(4,[(0,1),(1,2),(2,3)]),(4,[(0,1),(0,2),(0,3)]),(4,[(0,1),(2,3)])])
    patterns=checks=supports=0
    for n,edges in fixtures:
        p,c,s=fixture(n,edges,rng);patterns+=p;checks+=c;supports+=s
    print(f'PASS: {len(fixtures)} actual control graphs; {patterns} source-pattern algorithm outputs; {checks} direct state/input matching and reachability checks; {supports} robust assignment supports in exact-optimum enumerations.')
if __name__=='__main__':main()
