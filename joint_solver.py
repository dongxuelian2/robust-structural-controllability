"""Cubic-time 7/4 approximation by joint repair of a quota-spanning placement.

Input: nonempty iterable of outgoing adjacency bitmasks. Unit-cost, distinct
state dedicated actuators; return None precisely for infeasible instances.
See joint-repair-proof.md for the theorem, prior work and returned certificates.
"""
from collections import deque
from two_approx_solver import solve as baseline

def solve(rows):
    rows=tuple(rows)
    result=baseline(rows)
    if result is None: return None
    n=len(rows); base=set(result['relaxation_actuators'])
    critical=set(result['critical_states']); mate=dict(result['matching'])
    coverage={e:set() for e in range(n) if e not in base}
    incoming=[[i for i,row in enumerate(rows) if row>>j&1] for j in range(n)]
    for b in critical:
        q=deque([b]); reached={b}; lefts=set()
        while q:
            right=q.popleft()
            for left in incoming[right]:
                if left in lefts: continue
                lefts.add(left)
                assert left in mate
                e=mate[left]
                if e not in reached:
                    reached.add(e); q.append(e)
        for e in reached-{b}: coverage[e].add(b)
    uncovered=set(critical); repairs=[]
    while uncovered:
        e=max(coverage,key=lambda e:(len(coverage[e]&uncovered),-e))
        if len(coverage[e]&uncovered)<2: break
        repairs.append(e); uncovered-=coverage[e]
    residual=set(uncovered)
    for b in sorted(residual):
        e=min(e for e in coverage if b in coverage[e])
        repairs.append(e)
    result.update(actuators=tuple(sorted(base|set(repairs))),repair_states=tuple(repairs),residual=tuple(sorted(residual)),coverage={e:tuple(sorted(v)) for e,v in coverage.items()})
    return result
