from collections import Counter
from random import Random
import json
import two_approx_solver as a
import verify_core_theorem as v
v.CHECK_DIRECT=False
counts=Counter(); worst=(0,None)
def check(rows):
    global worst
    n=len(rows); full=(1<<n)-1
    result=a.solve(rows)
    d,opt=v.check(rows)
    counts['graphs']+=1
    if opt is None:
        assert result is None
        counts['infeasible']+=1; return
    assert result is not None
    counts['feasible']+=1
    source=v.sources(rows)
    relaxation=min(s.bit_count() for s in range(1<<n)
                   if all((s&c).bit_count()>=2 for c in source)
                   and v.matching(rows,full^s)==n-s.bit_count())
    assert result['relaxation_value']==relaxation
    assert len(result['critical_states'])<=d
    assert len(result['repair_states'])<=len(result['critical_states'])
    out=result['actuators']; mask=sum(1<<i for i in out)
    for i in out:
        remaining=mask^(1<<i)
        assert all(remaining&c for c in source)
        assert v.matching(rows,full^remaining)==n-remaining.bit_count()
    assert len(out)<=relaxation+d<=opt+d
    assert len(out)<=2*opt-(1 if d else 0)
    if d==0: assert len(out)==opt
    if len(out)/opt>worst[0]: worst=(len(out)/opt,{'rows':rows,'d':d,'opt':opt,'result':result})
for n in range(1,5):
    for bits in range(1<<(n*n)):
        check([(bits>>(n*i))&((1<<n)-1) for i in range(n)])
rng=Random(20261103)
for n in range(5,9):
    for _ in range(150):
        prob=rng.choice([.1,.2,.35,.5,.8])
        check([sum(1<<j for j in range(n) if rng.random()<prob) for i in range(n)])
print(json.dumps({'counts':dict(counts),'worst':worst,'result':'PASS'},indent=2))
