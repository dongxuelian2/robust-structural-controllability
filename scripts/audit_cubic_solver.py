"""Independent flow-relaxation and coloop-repair audit; no external packages."""
from collections import deque, Counter
from random import Random
from two_approx_solver import solve as cubic_solve
from audit_bounded_solver import matching_sets, oracle_source_sccs

def relaxation(rows,cs):
 n=len(rows); s=2*n+len(cs); t=s+1; N=t+1
 cap=[[0]*N for _ in range(N)]; adj=[set() for _ in range(N)]
 def edge(a,b,c):
  cap[a][b]=c; adj[a].add(b); adj[b].add(a)
 for u in range(n):
  edge(s,u,1)
  for v in range(n):
   if rows[u]>>v&1: edge(u,n+v,1)
 for v in range(n):
  group=next((i for i,c in enumerate(cs) if c>>v&1),None)
  edge(n+v,2*n+group if group is not None else t,1)
 for i,c in enumerate(cs): edge(2*n+i,t,c.bit_count()-2)
 flow=0
 while True:
  par=[-1]*N; par[s]=s; q=deque([s])
  while q and par[t]<0:
   u=q.popleft()
   for v in sorted(adj[u]):
    if cap[u][v]>0 and par[v]<0: par[v]=u; q.append(v)
  if par[t]<0: break
  v=t
  while v!=s:
   u=par[v];cap[u][v]-=1;cap[v][u]+=1;v=u
  flow+=1
 I=sum(1<<v for v in range(n) if any(cap[n+v][u] for u in range(n)))
 assert I.bit_count()==flow
 return ((1<<n)-1)^I

def check(rows):
 n=len(rows); full=(1<<n)-1; cs=oracle_source_sccs(rows)
 actual=cubic_solve(rows)
 matchable=matching_sets(rows); nu=max(s.bit_count() for s in matchable);d=n-nu
 if any(c.bit_count()<2 for c in cs):
  assert not any(full^s in matchable and all((s&c).bit_count()>=2 for c in cs) for s in range(full+1))
  assert actual is None
  return False,0,0,d
 ranks=[s.bit_count()+max((m.bit_count() for m in matchable if not m&s),default=0)-nu for s in range(full+1)]
 def spanning(s):return ranks[s]==d
 def robust(s):return spanning(s) and all(spanning(s^(1<<v)) for v in range(n) if s>>v&1)
 quota=lambda s:all((s&c).bit_count()>=2 for c in cs)
 opt=min(s.bit_count() for s in range(full+1) if quota(s) and robust(s))
 L=min(s.bit_count() for s in range(full+1) if quota(s) and spanning(s))
 assert actual is not None
 base=sum(1<<v for v in actual['relaxation_actuators'])
 chosen=sum(1<<v for v in actual['actuators'])
 assert base.bit_count()==L==actual['relaxation_value'] and spanning(base) and quota(base)
 assert chosen&base==base and robust(chosen) and quota(chosen)
 critical={v for v in range(n) if base>>v&1 and not spanning(base^(1<<v))}
 assert set(actual['critical_states'])==critical and len(critical)<=d
 repairs=sum(1<<v for v in actual['repair_states'])
 assert not repairs&base and chosen==base|repairs
 assert repairs.bit_count()<=len(critical)
 assert chosen.bit_count()<=L+d<=2*opt
 if d:assert chosen.bit_count()<=2*opt-1
 pairs=actual['matching']
 assert len({u for u,v in pairs})==len(pairs)==len({v for u,v in pairs})
 assert all(rows[u]>>v&1 for u,v in pairs)
 assert sum(1<<v for u,v in pairs)==full^base
 for b in critical:
  assert any(spanning((base^(1<<b))|1<<e) for e in actual['repair_states'])
 p=relaxation(rows,cs)
 assert p.bit_count()==L and quota(p) and spanning(p),(rows,'flow',p,L)
 initial=p;steps=0
 while True:
  coloops=[v for v in range(n) if p>>v&1 and not spanning(p^(1<<v))]
  assert len(coloops)<=d
  if not coloops:break
  b=coloops[0]
  e=next((v for v in range(n) if not p>>v&1 and spanning((p^(1<<b))|1<<v)),None)
  assert e is not None,(rows,'missing repair')
  old=p;p|=1<<e;steps+=1
  new=[v for v in range(n) if p>>v&1 and not spanning(p^(1<<v))]
  assert set(new)<set(coloops) and b not in new
 assert quota(p) and robust(p) and steps<=d
 assert p.bit_count()<=L+d<=2*opt
 if d:assert p.bit_count()<=2*opt-1
 else:assert p.bit_count()==opt
 return True,p.bit_count(),opt,d

if __name__=='__main__':
 count=Counter();worst=(0,None); examples=[]
 for n in range(1,5):
  for graph in range(1<<(n*n)):
   rows=tuple((graph>>(n*i))&((1<<n)-1) for i in range(n))
   f,a,o,d=check(rows);count['exhaustive_graphs']+=1;count['feasible']+=f
   if f and a/o>worst[0]:worst=(a/o,(rows,a,o,d))
  print('completed n=',n,flush=True)
 rng=Random(20261103)
 for n in range(5,9):
  for trial in range(100):
   density=rng.choice([.12,.25,.5,.75]);rows=tuple(sum(1<<v for v in range(n) if rng.random()<density) for u in range(n))
   f,a,o,d=check(rows);count['random_graphs']+=1;count['random_feasible']+=f
   if f and a/o>worst[0]:worst=(a/o,(rows,a,o,d))
 print(dict(count));print('largest observed ratio / instance:',worst)
 print('PASS: cubic solver exact relaxation, exact coloop set, matching certificates, repair witnesses, robust outputs, and size guarantee; plus independent flow/repair audit')
