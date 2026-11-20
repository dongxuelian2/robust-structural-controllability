from three_halves_solver import maximum_matching
from functools import lru_cache
from random import Random

def check(adj):
 n=len(adj);bits=[sum(1<<v for v in ns) for ns in adj]
 @lru_cache(None)
 def dp(mask):
  if not mask:return 0
  bit=mask&-mask;v=bit.bit_length()-1;rest=mask^bit
  return max([dp(rest)]+[1+dp(rest^(1<<u)) for u in adj[v] if rest>>u&1])
 before=[list(ns) for ns in adj]
 m=maximum_matching(adj)
 assert adj==before
 assert all(v==-1 or (v in adj[u] and m[v]==u) for u,v in enumerate(m))
 assert sum(v>=0 for v in m)==2*dp((1<<n)-1),(adj,m,dp((1<<n)-1))
 return len(m)
if __name__=='__main__':
 count=0
 for n in range(7):
  pairs=[(u,v) for u in range(n) for v in range(u+1,n)]
  for bits in range(1<<len(pairs)):
   adj=[[] for _ in range(n)]
   for i,(u,v) in enumerate(pairs):
    if bits>>i&1:adj[u].append(v);adj[v].append(u)
   check(adj);count+=1
  print('exhaustive',n,count,flush=True)
 rng=Random(202611040530)
 for n in range(7,19):
  for _ in range(100):
   adj=[[] for _ in range(n)];p=rng.choice([.1,.3,.5,.8])
   for u in range(n):
    for v in range(u+1,n):
     if rng.random()<p:adj[u].append(v);adj[v].append(u)
   check(adj);count+=1
  print('random',n,count,flush=True)
 print('PASS',count)
