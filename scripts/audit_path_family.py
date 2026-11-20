"""Independent construction and explicit witnesses for the 5/3 greedy family.
Run beside joint_solver.py: python3 audit_path_family.py
Only the candidate algorithm is imported; all family witnesses are explicit.
"""
from joint_solver import solve

def verify_matching(rows,pairs,right_set):
 assert len({u for u,v in pairs})==len(pairs)==len({v for u,v in pairs})
 assert all(rows[u]>>v&1 for u,v in pairs)
 assert {v for u,v in pairs}==right_set

for d in range(3,40,3):
 order=list(range(1,d-1,3))+[j for j in range(d-1) if j%3!=1]
 q={j:i for i,j in enumerate(order)}
 c=list(range(d-1,2*d-2));p=list(range(2*d-2,3*d-2));n=3*d-2
 rows=[0]*n
 for j in range(d-1): rows[c[j]]=sum(1<<v for v in (p[j],p[j+1],q[j]))
 for v in list(q.values())+p: rows[v]=sum(1<<u for u in c)
 # Strong connectivity, independently by graph search from each vertex.
 for start in range(n):
  reached={start};queue=[start]
  for u in queue:
   for v in range(n):
    if rows[u]>>v&1 and v not in reached:reached.add(v);queue.append(v)
  assert len(reached)==n
 F=[(c[j],q[j]) for j in range(d-1)]+[(q[j],c[j]) for j in range(d-1)]
 verify_matching(rows,F,set(range(n))-set(p))
 assert len(F)==2*d-2
 # Upper bound: each edge joins left centers to right noncenters or vice versa,
 # so each side contributes at most d-1. Thus deficiency is exactly d.
 assert all(not rows[u]&sum(1<<v for v in c) for u in c)
 assert all(not rows[u]&sum(1<<v for v in list(q.values())+p) for u in list(q.values())+p)
 W={p[0],p[-1]}|set(q.values());assert len(W)==d+1
 for deleted in W:
  targets=[]
  if deleted==p[0]:targets=[(c[j],p[j]) for j in range(d-1)]
  elif deleted==p[-1]:targets=[(c[j],p[j+1]) for j in range(d-1)]
  else:
   j=next(j for j in q if q[j]==deleted)
   targets=[(c[j],q[j])]+[(c[k],p[k+1]) for k in range(j)]+[(c[k],p[k]) for k in range(j+1,d-1)]
  pairs=targets+[(q[j],c[j]) for j in range(d-1)]
  verify_matching(rows,pairs,set(range(n))-(W-{deleted}))
 answer=solve(rows)
 assert set(answer['relaxation_actuators'])==set(p)
 assert len(answer['actuators'])==5*d//3 and len(answer['residual'])==d//3
 assert all(set(answer['coverage'][q[j]])=={p[j],p[j+1]} for j in q)
 print(f'd={d}: deficiency={d}, certified OPT={d+1}, actual output={5*d//3}')
print('PASS: explicit matching witnesses, deficiency certificate, robust optimum certificate, and actual greedy 5/3 family')
