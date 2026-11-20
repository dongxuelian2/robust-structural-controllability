"""Actual joint solver vs independent subset-DP ranks and optimal robust sets."""
from joint_solver import solve
from audit_seven_fourths import data
from collections import Counter
from random import Random

def check(rows):
 n=len(rows);full=(1<<n)-1;rank,d,cs=data(rows);answer=solve(iter(rows))
 if any(c.bit_count()<2 for c in cs):
  assert answer is None;return False
 assert answer is not None
 quota=lambda s:all((s&c).bit_count()>=2 for c in cs)
 robust=lambda s:rank[s]==d and all(rank[s^(1<<i)]==d for i in range(n) if s>>i&1)
 L=min(s.bit_count() for s in range(full+1) if rank[s]==d and quota(s))
 opt=min(s.bit_count() for s in range(full+1) if robust(s) and quota(s))
 P=sum(1<<i for i in answer['relaxation_actuators']);S=sum(1<<i for i in answer['actuators'])
 assert P.bit_count()==L==answer['relaxation_value'] and rank[P]==d and quota(P)
 D={i for i in range(n) if P>>i&1 and rank[P^(1<<i)]<d}
 assert set(answer['critical_states'])==D
 expected={e:tuple(sorted(b for b in D if rank[(P^(1<<b))|1<<e]==d)) for e in range(n) if not P>>e&1}
 assert answer['coverage']==expected,(rows,'coverage',expected,answer['coverage'])
 residual=set(answer['residual']);r=len(residual)
 cuts=[{b}|{e for e,A in expected.items() if b in A} for b in residual]
 assert all(not cuts[i]&cuts[j] for i in range(r) for j in range(i))
 assert opt>=2*r
 repairs=answer['repair_states'];assert len(repairs)==len(set(repairs))
 assert 2*len(repairs)<=len(D)+r
 assert robust(S) and quota(S) and 4*S.bit_count()<=7*opt
 return True

counts=Counter()
for n in range(1,5):
 for graph in range(1<<(n*n)):
  rows=tuple((graph>>(n*i))&((1<<n)-1) for i in range(n))
  counts['exhaustive_graphs']+=1;counts['exhaustive_feasible']+=check(rows)
 print('completed n=',n,flush=True)
rng=Random(20261104)
for n in range(5,9):
 for trial in range(200):
  p=rng.choice([.12,.25,.5,.75]);rows=tuple(sum(1<<j for j in range(n) if rng.random()<p) for i in range(n))
  counts['random_graphs']+=1;counts['random_feasible']+=check(rows)
print(dict(counts));print('PASS: actual implementation, full repair sets, exact relaxation, residual disjointness, robust output, 7/4, and iterator inputs')
