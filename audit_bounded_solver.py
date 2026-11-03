"""Independent audit: compare returned placements with exact subset-DP oracle."""
from collections import Counter
from random import Random
from bounded_core_solver import solve, matching_number, source_components
def matching_sets(out):
    """All matchable right subsets, via exact left-vertex dynamic programming."""
    masks = {0}
    for neighbors in out:
        for mask in tuple(masks):
            available = neighbors & ~mask
            while available:
                bit = available & -available
                masks.add(mask | bit)
                available -= bit
    return masks



def oracle_source_sccs(out):
    """Independent source-SCC check: per-vertex DFS reachability."""
    n = len(out)
    reach = []
    for start in range(n):
        seen, stack = {start}, [start]
        while stack:
            u = stack.pop()
            for v in range(n):
                if out[u] >> v & 1 and v not in seen:
                    seen.add(v)
                    stack.append(v)
        reach.append(seen)
    components = {frozenset(v for v in range(n)
                           if v in reach[u] and u in reach[v]) for u in range(n)}
    return sorted(sum(1 << v for v in c) for c in components
                  if not any(out[u] >> v & 1 for u in range(n) if u not in c for v in c))




def audit(rows):
    n=len(rows); full=(1<<n)-1
    matchable=matching_sets(rows); cs=oracle_source_sccs(rows)
    assert sorted(source_components(rows))==cs
    nu=max(s.bit_count() for s in matchable)
    assert matching_number(rows,full)==nu
    ordinary=[full^s in matchable and all(s&c for c in cs) for s in range(full+1)]
    robust=[ordinary[s] and all(ordinary[s^(1<<i)] for i in range(n) if s>>i&1) for s in range(full+1)]
    optimum=min((s.bit_count() for s in range(full+1) if robust[s]),default=None)
    before=tuple(rows); answer=solve(rows)
    assert tuple(rows)==before
    if answer is None:
        assert optimum is None,(rows,answer,optimum)
    else:
        assert len(set(answer))==len(answer)
        selected=sum(1<<i for i in answer)
        assert robust[selected] and len(answer)==optimum,(rows,answer,optimum)
    return n-nu, optimum is not None

if __name__=='__main__':
    counts=Counter(); by_d=Counter()
    for n in range(1,5):
        for bits in range(1<<(n*n)):
            rows=tuple((bits>>(n*i))&((1<<n)-1) for i in range(n))
            d,feasible=audit(rows)
            counts['exhaustive_graphs']+=1; counts['exhaustive_feasible']+=feasible; by_d[d]+=1
        print('completed n=',n,flush=True)
    rng=Random(20261103)
    for n in range(5,9):
        for trial in range(100):
            p=rng.choice([.12,.25,.5,.75])
            rows=tuple(sum(1<<j for j in range(n) if rng.random()<p) for i in range(n))
            d,feasible=audit(rows); counts['random_graphs']+=1; counts['random_feasible']+=feasible
    for rows in ((),(-1,), (2,), ('x',)):
        try: solve(rows)
        except ValueError: pass
        else: raise AssertionError(('invalid input accepted',rows))
    print(dict(counts)); print('exhaustive deficiency histogram:',dict(sorted(by_d.items())))
    print('PASS: solver placements, optimality, infeasibility, nonmutation, and invalid inputs')
