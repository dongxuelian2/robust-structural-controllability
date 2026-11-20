"""Independent exhaustive oracle audit of fixed-P greedy repair, 7/4 bound."""
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

def check(rows,all_optimal=False):
 n=len(rows);full=(1<<n)-1;rank,d,cs=data(rows)
 if any(c.bit_count()<2 for c in cs):return None
 quota=lambda s: all((s&c).bit_count()>=2 for c in cs)
 robust=lambda s: rank[s]==d and all(rank[s^(1<<i)]==d for i in range(n) if s>>i&1)
 opt=min(s.bit_count() for s in range(full+1) if quota(s) and robust(s))
 possible=[s for s in range(full+1) if quota(s) and rank[s]==d]
 L=min(s.bit_count() for s in possible)
 ps=[s for s in possible if s.bit_count()==L]
 worst=(0,None)
 for P in (ps if all_optimal else ps[:1]):
  D={b for b in range(n) if P>>b&1 and rank[P^(1<<b)]<d}
  assert len(D)<=d
  A={e:{b for b in D if rank[(P^(1<<b))|1<<e]==d} for e in range(n) if not P>>e&1}
  remaining=set(D);chosen=[];greedy=0
  while True:
   e=next((e for e,a in A.items() if len(a&remaining)>=2),None)
   if e is None:break
   chosen.append(e);greedy+=1;remaining-=A[e]
  r=len(remaining)
  cuts=[{b}|{e for e,a in A.items() if b in a} for b in remaining]
  assert all(not cuts[i]&cuts[j] for i in range(r) for j in range(i))
  # Each of these sets is the complement of a rank-(d-1) hyperplane.
  for b,cut in zip(remaining,cuts):
   mask=sum(1<<i for i in cut);H=full^mask
   assert rank[H]==d-1 and all(rank[H|1<<i]==d for i in cut)
  assert opt>=2*r
  for b in remaining:
   e=next(e for e,a in A.items() if b in a)
   chosen.append(e)
  assert len(set(chosen))==len(chosen)
  S=P|sum(1<<e for e in chosen)
  assert robust(S) and quota(S)
  assert 2*len(chosen)<=len(D)+r
  assert 4*S.bit_count()<=7*opt
  if d:assert 4*S.bit_count()<=7*opt-2
  else:assert S.bit_count()==opt
  candidate=(S.bit_count()/opt,(rows,P,S,opt,d,len(D),r))
  if candidate[0]>worst[0]:worst=candidate
 return worst

if __name__=='__main__':
 counts=Counter();worst=(0,None)
 for n in range(1,5):
  for bits in range(1<<(n*n)):
   rows=tuple((bits>>(n*i))&((1<<n)-1) for i in range(n))
   result=check(rows,all_optimal=True);counts['exhaustive_graphs']+=1
   if result:
    counts['exhaustive_feasible']+=1
    if result[0]>worst[0]:worst=result
  print('completed n=',n,flush=True)
 rng=Random(20261104)
 for n in range(5,9):
  for _ in range(200):
   p=rng.choice([.12,.25,.5,.75])
   rows=tuple(sum(1<<j for j in range(n) if rng.random()<p) for i in range(n))
   result=check(rows,all_optimal=True);counts['random_graphs']+=1
   if result:
    counts['random_feasible']+=1
    if result[0]>worst[0]:worst=result
 print(dict(counts));print('worst observed:',worst)
 print('PASS: all optimum relaxation sets, greedy repair, disjoint cocircuits, residual lower bound, robustness, and 7/4 guarantee')
