"""Independent exhaustive subset-DP feasibility/optimality oracle."""
from three_halves_solver import solve
from audit_bounded_solver import matching_sets,oracle_source_sccs
from functools import lru_cache
from random import Random
import json


def matching(covers, D):
 """Exact exponential matching oracle, independent of blossom implementation."""
 neighbors={b:0 for b in range(D.bit_length()) if D>>b&1}
 for covered in covers.values():
  for b in neighbors:
   if covered>>b&1:neighbors[b]|=covered^(1<<b)
 @lru_cache(None)
 def optimum(U):
  if not U:return 0
  v=(U&-U).bit_length()-1;rest=U^(1<<v)
  return max([optimum(rest)]+[1+optimum(rest^(1<<w))
             for w in neighbors if rest>>w&1 and neighbors[v]>>w&1])
 return optimum(D)

count=0;feasible=0;worst=(0,None)
def check(rows):
 global count,feasible,worst
 count+=1;before=tuple(rows);result=solve(rows)
 assert tuple(rows)==before
 n=len(rows);full=(1<<n)-1;cs=oracle_source_sccs(rows)
 ms=matching_sets(rows)
 span=[full^s in ms for s in range(full+1)]
 quota=[all((s&c).bit_count()>=2 for c in cs) for s in range(full+1)]
 robust=[span[s] and quota[s] and all(span[s^(1<<i)] for i in range(n) if s>>i&1) for s in range(full+1)]
 opt=min((s.bit_count() for s in range(full+1) if robust[s]),default=None)
 if opt is None:
  assert result is None
  return
 feasible+=1;assert result is not None
 L=min(s.bit_count() for s in range(full+1) if quota[s] and span[s])
 P=sum(1<<b for b in result['relaxation_actuators'])
 S=sum(1<<b for b in result['actuators'])
 D=sum(1<<b for b in result['critical_states'])
 assert P.bit_count()==L==result['relaxation_value'] and quota[P] and span[P]
 assert D==sum(1<<b for b in range(n) if P>>b&1 and not span[P^(1<<b)])
 covers={e:sum(1<<b for b in range(n) if D>>b&1 and span[(P^(1<<b))|1<<e]) for e in range(n) if not P>>e&1}
 assert covers=={e:sum(1<<b for b in bs) for e,bs in result['repair_supports']}
 mu=matching(covers,D)
 assert mu==result['pair_matching_size']
 used=set()
 for b,c,e in result['pair_repairs']:
  assert b not in used and c not in used and b!=c
  used.update((b,c));assert covers[e]>>b&1 and covers[e]>>c&1
 for b,e in result['singleton_repairs']:
  assert b not in used;used.add(b);assert covers[e]>>b&1
 assert used==set(result['critical_states'])
 assert S&P==P and robust[S]
 assert S.bit_count()==len(set(result['actuators']))
 q=D.bit_count()-mu
 assert q==result['repair_budget'] and 2*q<=opt
 assert result['opt_lower_bound']<=opt
 assert S.bit_count()<=L+q and 2*S.bit_count()<=3*opt
 pairs=result['matching']
 assert len({u for u,v in pairs})==len({v for u,v in pairs})==len(pairs)
 assert all(rows[u]>>v&1 for u,v in pairs)
 assert sum(1<<v for u,v in pairs)==full^P
 if S.bit_count()/opt>worst[0]:worst=(S.bit_count()/opt,dict(rows=rows,opt=opt,result=result))

if __name__=='__main__':
 for n in range(1,5):
  for bits in range(1<<(n*n)):
   check(tuple(bits>>(n*i)&((1<<n)-1) for i in range(n)))
  print('exhaustive',n,count,feasible,flush=True)
 rng=Random(202611040531)
 for n in range(5,12):
  for _ in range(500):
   p=rng.choice([.1,.15,.2,.3,.5,.8]);check(tuple(sum(1<<j for j in range(n) if rng.random()<p) for i in range(n)))
  print('random',n,count,feasible,flush=True)
 for rows in [(),(-1,), (2,),('x',)]:
  try:solve(rows)
  except ValueError:pass
  else:raise AssertionError(('invalid input accepted',rows))
 print('PASS',count,feasible,json.dumps(worst),flush=True)
