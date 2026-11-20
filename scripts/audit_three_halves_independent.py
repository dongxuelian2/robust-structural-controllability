"""Independent exact-oracle audit of co-cover matching 3/2 approximation.
Uses subset-DP transversal ranks and subset-DP general graph matching.
The matching oracle here is exponential and is only a small-instance verifier.
"""
from three_halves_solver import solve as actual_solve
from functools import lru_cache
from collections import Counter
from random import Random

def data(rows):
 n=len(rows); full=(1<<n)-1
 matchable={0}
 for row in rows:
  for s in tuple(matchable):
   avail=row&~s
   while avail:
    bit=avail&-avail;avail-=bit;matchable.add(s|bit)
 rt=[s.bit_count() if s in matchable else 0 for s in range(full+1)]
 for i in range(n):
  for s in range(full+1):
   if s>>i&1: rt[s]=max(rt[s],rt[s^(1<<i)])
 nu=rt[full];d=n-nu
 rank=[s.bit_count()+rt[full^s]-nu for s in range(full+1)]
 reach=[row|1<<i for i,row in enumerate(rows)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 remaining=full;cs=[]
 while remaining:
  i=(remaining&-remaining).bit_length()-1
  c=sum(1<<j for j in range(n) if reach[i]>>j&1 and reach[j]>>i&1)
  remaining&=~c
  if not any(row&c for j,row in enumerate(rows) if not c>>j&1):cs.append(c)
 return rank,d,cs


def check(rows):
 n=len(rows);full=(1<<n)-1;rank,d,cs=data(rows)
 actual=actual_solve(iter(rows))
 if any(c.bit_count()<2 for c in cs):
  assert actual is None
  return None
 quota=lambda s:all((s&c).bit_count()>=2 for c in cs)
 robust=lambda s:rank[s]==d and all(rank[s^(1<<i)]==d for i in range(n) if s>>i&1)
 opt=min(s.bit_count() for s in range(full+1) if quota(s) and robust(s))
 candidates=[s for s in range(full+1) if quota(s) and rank[s]==d]
 L=min(s.bit_count() for s in candidates);worst=(0,None)
 assert actual is not None
 actualP=sum(1<<b for b in actual['relaxation_actuators'])
 actualS=sum(1<<b for b in actual['actuators'])
 assert actualP.bit_count()==L==actual['relaxation_value'] and quota(actualP) and rank[actualP]==d
 assert actualS&actualP==actualP and quota(actualS) and robust(actualS)
 assert 2*actualS.bit_count()<=3*opt
 assert actual['opt_lower_bound']<=opt
 for P in (s for s in candidates if s.bit_count()==L):
  D=[b for b in range(n) if P>>b&1 and rank[P^(1<<b)]<d];k=len(D)
  A={e:sum(1<<j for j,b in enumerate(D) if rank[(P^(1<<b))|1<<e]==d) for e in range(n) if not P>>e&1}
  cuts=[(1<<b)|sum(1<<e for e,a in A.items() if a>>j&1) for j,b in enumerate(D)]
  adj=[0]*k
  for a in A.values():
   for j in range(k):
    if a>>j&1:adj[j]|=a^(1<<j)
  @lru_cache(None)
  def maximum(s):
   if not s:return 0,()
   j=(s&-s).bit_length()-1;rest=s^(1<<j);best=maximum(rest)
   neighbors=adj[j]&rest
   while neighbors:
    bit=neighbors&-neighbors;neighbors-=bit;i=bit.bit_length()-1
    val,pairs=maximum(rest^bit)
    if val+1>best[0]:best=(val+1,pairs+((j,i),))
   return best
  nu,pairs=maximum((1<<k)-1)
  if P==actualP:
   assert tuple(D)==actual['critical_states']
   supports={e:tuple(D[j] for j in range(k) if a>>j&1) for e,a in A.items()}
   assert dict(actual['repair_supports'])==supports
   assert actual['pair_matching_size']==nu and actual['repair_budget']==k-nu
   assert len(actual['repair_states'])<=k-nu
   covered=set()
   for b,c,e in actual['pair_repairs']:
    assert b!=c and b not in covered and c not in covered
    covered|={b,c};assert b in supports[e] and c in supports[e]
   for b,e in actual['singleton_repairs']:
    assert b not in covered and b in supports[e];covered.add(b)
   assert covered==set(D)
   cert=actual['matching']
   assert len({u for u,v in cert})==len(cert)==len({v for u,v in cert})
   assert all(rows[u]>>v&1 for u,v in cert)
   assert {v for u,v in cert}=={v for v in range(n) if not P>>v&1}
  # Check the critical union-rank identity for EVERY subset of original D.
  for R in range(1<<k):
   union=0
   for j in range(k):
    if R>>j&1:union|=cuts[j]
   assert rank[full^union]==d-R.bit_count(),(rows,P,R,'rank identity')
  tb=-10**9
  for U in range(1<<k):
   unseen=((1<<k)-1)^U;parts=[]
   while unseen:
    bit=unseen&-unseen;component=bit;front=bit;unseen^=bit
    while front:
     bit=front&-front;front-=bit;j=bit.bit_length()-1
     add=adj[j]&unseen;component|=add;front|=add;unseen^=add
    parts.append(component)
   unions=[]
   for component in parts:
    union=0
    for j in range(k):
     if component>>j&1:union|=cuts[j]
    unions.append(union)
   assert all(not unions[i]&unions[j] for i in range(len(unions)) for j in range(i))
   assert opt>=k-U.bit_count()+len(parts)
   odd=sum(s.bit_count()%2 for s in parts)
   tb=max(tb,odd-U.bit_count())
  assert 2*nu==k-tb
  assert opt>=2*(k-nu),(rows,P,opt,k,nu,'packing bound')
  repairs=set();matched=0
  for i,j in pairs:
   repairs.add(next(e for e,a in A.items() if a>>i&1 and a>>j&1));matched|=(1<<i)|(1<<j)
  for j in range(k):
   if not matched>>j&1:repairs.add(next(e for e,a in A.items() if a>>j&1))
  S=P|sum(1<<e for e in repairs)
  assert len(repairs)<=k-nu
  assert quota(S) and robust(S) and 2*S.bit_count()<=3*opt
  observed=(S.bit_count()/opt,(rows,P,S,opt,d,k,nu))
  if observed[0]>worst[0]:worst=observed
 return worst

if __name__=='__main__':
 counts=Counter();worst=(0,None)
 for n in range(1,5):
  for graph in range(1<<(n*n)):
   rows=tuple((graph>>(n*i))&((1<<n)-1) for i in range(n))
   result=check(rows);counts['exhaustive_graphs']+=1
   if result:
    counts['exhaustive_feasible']+=1
    if result[0]>worst[0]:worst=result
  print('completed n=',n,flush=True)
 rng=Random(20261104)
 for n in range(5,9):
  for trial in range(200):
   p=rng.choice([.12,.25,.5,.75]);rows=tuple(sum(1<<j for j in range(n) if rng.random()<p) for i in range(n))
   result=check(rows);counts['random_graphs']+=1
   if result:
    counts['random_feasible']+=1
    if result[0]>worst[0]:worst=result
 print(dict(counts));print('worst observed:',worst)
 print('PASS: actual blossom solver vs independent full ranks/general-matching DP; supports, exact pairs, certificates, iterable input, robustness and 3/2; plus all-P mathematical checks')
