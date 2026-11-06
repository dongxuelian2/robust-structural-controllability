"""Independent full-source-copy, path-cover, and approximation-transfer audit.

Uses resource/job matchings, not the bicircular degree criterion.
"""
import itertools

def isforest(n,es):
 p=list(range(n));d=[0]*n
 def find(v):
  while p[v]!=v:v=p[v]
  return v
 for a,b in es:
  d[a]+=1;d[b]+=1
  if d[a]>2 or d[b]>2 or find(a)==find(b):return False
  p[find(a)]=find(b)
 return True

def matching(n,res,omit):
 at=[[] for _ in range(n)]
 for j,r in enumerate(res):
  if j!=omit:
   for v in r:at[v].append(j)
 mate={}
 def aug(v,seen):
  for j in at[v]:
   if j in seen:continue
   seen.add(j)
   if j not in mate or aug(mate[j],seen):mate[j]=v;return True
  return False
 return all(aug(v,set()) for v in range(n))
def robust(n,res):
 return len(res)>=n+1 and all(matching(n,res,j) for j in range(len(res)))
def degs(n,es):
 d=[0]*n
 for a,b in es:d[a]+=1;d[b]+=1
 return d

def anchored_pc(n,es,s,t):
 best=n
 for mask in range(1<<len(es)):
  fs=[e for j,e in enumerate(es) if mask>>j&1]
  d=degs(n,fs)
  if d[s]<=1 and d[t]<=1 and isforest(n,fs):best=min(best,n-len(fs))
 return best

def fullopt(n,es,s,t):
 # completion singleton copies, physical endpoint copies, then source A/B
 rs=[set(e) for e in es]+[{v} for v in range(n)]+[{s},{t},{s},{t}]
 a,b=len(rs)-2,len(rs)-1
 for total in range(n+1,2*n+3):
  for k in range(max(n+1,total-2),min(total,len(rs))+1):
   for inds in itertools.combinations(range(len(rs)),k):
    fee=max(0,2-int(a in inds)-int(b in inds))
    if k+fee!=total:continue
    counts['supports']+=1
    if robust(n,[rs[i] for i in inds]):return total
 raise AssertionError('infeasible')

if __name__ == '__main__':
    counts={'fixtures':0,'supports':0,'covers':0,'good_ratio_covers':0}
    for n in range(2,5):
     pairs=list(itertools.combinations(range(n),2))
     for m in range(1<<len(pairs)):
      es=[e for j,e in enumerate(pairs) if m>>j&1]
      if any(all(v not in e for e in es) for v in range(n)):continue
      for s,t in itertools.combinations(range(n),2):
       p=anchored_pc(n,es,s,t)
       actual=fullopt(n,es,s,t)
       assert actual==n+p,(n,es,s,t,p,actual)
       ext=es+[(s,n),(t,n+1)]
       N=n+2;maxL=0
       for mask in range(1<<len(ext)):
        fs=[e for j,e in enumerate(ext) if mask>>j&1]
        if not isforest(N,fs):continue
        counts['covers']+=1
        L=len(fs);maxL=max(maxL,L)
        fs=set(fs)
        for v,w in ((s,n),(t,n+1)):
         if (v,w) in fs:continue
         inc=sorted(e for e in fs if v in e)
         if len(inc)==2:fs.remove(inc[0])
         fs.add((v,w))
        assert len(fs)>=L and isforest(N,fs)
        recovered=fs-{(s,n),(t,n+1)}
        dd=degs(n,recovered)
        q=n-len(recovered)
        assert dd[s]<=1 and dd[t]<=1 and q==N-len(fs)
        if 7*L<6*(N-p):continue
        counts['good_ratio_covers']+=1
        assert 7*(n+q)<=8*(n+p)
        # Lift edges+endpoint loops, and singleton+arbitrary original neighbor.
        resources=set(recovered)
        for v,d in enumerate(dd):
         if d<=1:resources.add((v,))
         if d==0:resources.add(next(e for e in es if v in e))
        assert (s,) in resources and (t,) in resources
        assert len(resources)<=n+q and robust(n,[set(e) for e in resources])
       assert N-maxL==p
       counts['fixtures']+=1
    print('INDEPENDENT ACTUAL-MATCHING ANCHORED AUDIT: ALL PASS',counts)
