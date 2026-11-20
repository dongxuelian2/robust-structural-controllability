#!/usr/bin/env python3
"""Finite audit of weighted cardinality quotas using the standard library."""
from itertools import combinations
from math import comb
from random import Random
import json

def subsets(E,b):
    for j in range(min(b,len(E))+1): yield from combinations(E,j)
def gf2rank(vs):
    piv={}
    for x in vs:
        while x:
            p=x.bit_length()-1
            if p in piv:x^=piv[p]
            else:piv[p]=x;break
    return len(piv)
def completion(C,w,col,beta):
    D=set(C)
    for g,q in beta.items():
        missing=max(0,q-sum(col[e]==g for e in D))
        D.update(sorted((e for e in range(len(w)) if col[e]==g and e not in D),key=lambda e:(w[e],e))[:missing])
    return sum(w[e] for e in D)
def adjusted(w,col,beta):
    a={g:sorted(w[e] for e in range(len(w)) if col[e]==g)[q-1] for g,q in beta.items()}
    return [max(0,w[e]-a[col[e]]) if col[e] in a else w[e] for e in range(len(w))]
def core(rank, r, k, w, col, safe, beta, *, prices=None, protected=()):
    b, H = (k + 1) * r, comb(r + k - 1, k)
    L, A = [0], [0]
    for _ in range(r):
        L.append(H * L[-1] + b)
        A.append(H * A[-1] + b * L[-1] + b)
    prices = adjusted(w, col, beta) if prices is None else prices
    output = set()
    discarded = {'color': 0, 'total': 0}
    regular_pool = [e for e in range(len(w)) if col[e] not in protected]
    protected_pools = [[e for e in range(len(w)) if col[e] == group] for group in protected]

    for pool_index, pool in enumerate([regular_pool, *protected_pools]):
        for status in (False, True):
            retained = []
            ordered = sorted((e for e in pool if safe[e] == status and rank([e]) > 0),
                             key=lambda e: (prices[e], w[e], e))
            for x in ordered:
                witness = None
                for basis in subsets(retained, r):
                    h = len(basis)
                    if rank(basis) != h or rank((*basis, x)) != h:
                        continue
                    flat = [e for e in retained if rank((*basis, e)) == h]
                    if sum(col[e] == col[x] for e in flat) == L[h]:
                        witness = 'color'
                        break
                    if pool_index == 0 and len(flat) == A[h]:
                        witness = 'total'
                        break
                if witness:
                    discarded[witness] += 1
                else:
                    retained.append(x)
            output.update(retained)
    return output, discarded

def resilient(C,rank,r,k,safe):
    V=[e for e in C if not safe[e]]
    return all(rank(set(C)-set(F))==r for F in subsets(V,k))

def main():
    rng = Random(16333)
    count={'marginal':0,'same_group':0,'optima':0,'exchange':0,'color_discards':0,'total_discards':0}
    # Independent marginal checks with direct completion enumeration formula.
    for z in range(150):
        n=rng.randint(3,10);w=[rng.randrange(7) for _ in range(n)];col=[rng.randrange(4) for _ in range(n)]
        beta={g:rng.randint(1,col.count(g)) for g in set(col) if g!=0};v=adjusted(w,col,beta)
        for Yt in subsets(list(range(n)),min(4,n)):
            Y=set(Yt);base=completion(Y,w,col,beta)
            for x in set(range(n))-Y:
                dx=completion(Y|{x},w,col,beta)-base
                assert dx>=v[x];count['marginal']+=1
                if col[x]==0 or not any(col[e]==col[x] for e in Y):assert dx==v[x]
                for y in set(range(n))-Y-{x}:
                    if col[y]==col[x] and w[y]<=w[x]:
                        assert completion(Y|{y},w,col,beta)<=completion(Y|{x},w,col,beta);count['same_group']+=1

    # General represented matroid small exhaustive optimum checks, loops retained for quotas.
    for z in range(80):
        n=rng.randint(4,9);d=rng.randint(1,3);vec=[rng.randrange(1<<d) for _ in range(n)]
        rank=lambda C:gf2rank(vec[e] for e in C)
        r=rank(range(n))
        if not r:continue
        k=rng.randint(0,2);w=[rng.randrange(6) for _ in range(n)];col=[rng.randrange(4) for _ in range(n)]
        BETA={g:rng.randint(1,col.count(g)) for g in set(col) if g};safe=[rng.choice([False,False,True]) for _ in range(n)]
        W,ds=core(rank,r,k,w,col,safe,BETA)
        full=min((sum(w[e] for e in B) for B in subsets(list(range(n)),n) if all(sum(col[e]==g for e in B)>=q for g,q in BETA.items()) and resilient(B,rank,r,k,safe)),default=None)
        reduced=min((completion(C,w,col,BETA) for C in subsets(sorted(W),(k+1)*r) if resilient(C,rank,r,k,safe)),default=None)
        assert full==reduced,(z,full,reduced);count['optima']+=1
        count['color_discards']+=ds['color'];count['total_discards']+=ds['total']

    # Rank-one instances force actual global saturation at low thresholds.
    for z in range(100):
        n=18;r=1;k=1;vec=[rng.choice([0,1,1,1]) for _ in range(n)]
        rank=lambda C:int(any(vec[e] for e in C));w=[rng.randrange(10) for _ in range(n)]
        col=[e//2+1 for e in range(n)];BETA={g:rng.randint(1,2) for g in set(col)};safe=[False]*n
        W,ds=core(rank,r,k,w,col,safe,BETA)
        full=min((completion(C,w,col,BETA) for C in subsets(list(range(n)),2) if resilient(C,rank,r,k,safe)),default=None)
        reduced=min((completion(C,w,col,BETA) for C in subsets(sorted(W),2) if resilient(C,rank,r,k,safe)),default=None)
        assert full==reduced;count['optima']+=1
        count['color_discards']+=ds['color'];count['total_discards']+=ds['total']

    # U_2,n with many colors triggers rank-two global saturation; audit all single replacements for sampled witnesses.
    n=105;r=2;k=1;rank=lambda C:min(2,len(set(C)));w=[rng.randrange(25) for _ in range(n)];col=[e//3+1 for e in range(n)]
    BETA={g:rng.randint(1,3) for g in set(col)};safe=[False]*n
    W,ds=core(rank,r,k,w,col,safe,BETA)
    count['color_discards']+=ds['color'];count['total_discards']+=ds['total']
    for z in range(2000):
        C=set(rng.sample(range(n),rng.choice([3,4])))
        for x in C-W:
            assert any(resilient(C-{x}|{y},rank,r,k,safe) and completion(C-{x}|{y},w,col,BETA)<=completion(C,w,col,BETA) for y in W-C)
            count['exchange']+=1
    print(json.dumps({'status':'passed','checks':count,'rank2_core_size':len(W),'rank2_input_size':n},indent=2))


if __name__ == '__main__':
    main()
