"""Independent finite audit. No representation assumption; uses bitmask ranks.

Run: python audit_weighted_core.py
The exhaustive checker enumerates all ground subsets and all allowed failures,
independently of the implementation's core, cardinality bound, or max-failure
shortcut. Weight enumeration is exhaustive only in the suites labeled so.
"""
from itertools import combinations, product
from functools import lru_cache
import json
import random
import time
from weighted_fault_core import solve, weighted_core

COUNTS = dict(cases=0, feasible=0, infeasible=0, capped_cases=0,
              rejected_cases=0, cap_and_rejection_cases=0, strict_core_cases=0,
              iterable_cases=0, checked_core_property=0,
              capped_rank2_k1_n_ge9=0, capped_rank2_k1_n_ge9_and_strict=0,
              capped_mixed_cases=0,
              capped_all_safe_cases=0, capped_all_vulnerable_cases=0)
SUITES = {}
EXAMPLES = {}


def submasks(m):
    s = m
    while True:
        yield s
        if s == 0:
            return
        s = (s - 1) & m


def make_rank(kind, n, **kwargs):
    full = (1 << n) - 1
    if kind == 'classes':
        labels = kwargs['labels']
        trunc = kwargs['trunc']
        def raw(m):
            return min(trunc, len({labels[i] for i in range(n)
                                  if m >> i & 1 and labels[i] is not None}))
    elif kind == 'sparse_paving':
        r = kwargs['r']
        forbidden = frozenset(sum(1 << i for i in s) for s in kwargs['forbidden'])
        assert all((a & b).bit_count() <= r-2 for a,b in combinations(forbidden,2))
        def raw(m):
            return r-1 if m in forbidden else min(r,m.bit_count())
    elif kind == 'graphic':
        edges = kwargs['edges']
        vertices = kwargs['vertices']
        def raw(m):
            parents = list(range(vertices))
            def root(v):
                while parents[v] != v:
                    v = parents[v]
                return v
            result = 0
            for i,(u,v) in enumerate(edges):
                if m >> i & 1:
                    a,b = root(u),root(v)
                    if a != b:
                        parents[a] = b
                        result += 1
            return result
    else:
        raise ValueError(kind)
    @lru_cache(None)
    def rankmask(m):
        return raw(m)
    def oracle(a):
        assert isinstance(a,frozenset), type(a)
        return rankmask(sum(1 << i for i in a))
    return rankmask, oracle


@lru_cache(None)
def members(n,m):
    return tuple(i for i in range(n) if m >> i & 1)


def feasible_masks(n, rankmask, k, vulnerable):
    r = rankmask((1 << n)-1)
    out = []
    for b in range(1 << n):
        if all(rankmask(b ^ failed) == r
               for failed in submasks(b & vulnerable) if failed.bit_count() <= k):
            out.append(b)
    return out


def run_case(name, weights, oracle, rankmask, k, vulnerable, feas=None,
             iterables=False, core_property=False):
    n = len(weights)
    if feas is None:
        feas = feasible_masks(n,rankmask,k,vulnerable)
    exact = min((sum(weights[i] for i in members(n,b)) for b in feas), default=None)
    v = members(n,vulnerable)
    args = ((w for w in weights),(x for x in v)) if iterables else (weights,v)
    result = solve(args[0], oracle, k, args[1])
    cost, witness, core, stats = result
    assert cost == exact, (name,weights,k,vulnerable,result,exact)
    if cost is not None:
        mask = sum(1 << i for i in witness)
        assert mask in feas, (name,'bad witness',witness)
        assert sum(weights[i] for i in witness) == cost
        assert witness <= core
    q = (k+1)*rankmask((1 << n)-1)
    assert witness is None or len(witness) <= q
    if core_property:
        cm = sum(1 << i for i in core)
        bysize = {}
        for b in feas:
            if b & ~cm == 0:
                size = b.bit_count()
                val = sum(weights[i] for i in members(n,b))
                bysize[size] = min(bysize.get(size,float('inf')),val)
        nonloops = sum(1 << i for i in range(n) if rankmask(1 << i) > 0)
        for b in feas:
            if b.bit_count() <= q and b & ~nonloops == 0:
                val = sum(weights[i] for i in members(n,b))
                assert bysize.get(b.bit_count(),float('inf')) <= val, (name,'strong cardinality',b)
        COUNTS['checked_core_property'] += 1
    COUNTS['cases'] += 1
    COUNTS['feasible' if exact is not None else 'infeasible'] += 1
    COUNTS['capped_cases'] += stats['capped'] > 0
    COUNTS['rejected_cases'] += stats['rejected'] > 0
    COUNTS['cap_and_rejection_cases'] += stats['capped'] > 0 and stats['rejected'] > 0
    COUNTS['strict_core_cases'] += len(core) < n
    COUNTS['iterable_cases'] += iterables
    COUNTS['capped_rank2_k1_n_ge9'] += (stats['capped'] > 0 and n >= 9 and k == 1
                                             and rankmask((1 << n)-1) == 2)
    COUNTS['capped_rank2_k1_n_ge9_and_strict'] += (stats['capped'] > 0 and n >= 9
                                                    and k == 1 and len(core) < n
                                                    and rankmask((1 << n)-1) == 2)
    if stats['capped']:
        status_key = ('capped_all_safe_cases' if vulnerable == 0 else
                      'capped_all_vulnerable_cases' if vulnerable == (1 << n)-1 else
                      'capped_mixed_cases')
        COUNTS[status_key] += 1
    SUITES[name] = SUITES.get(name,0)+1
    if name.startswith('adversarial'):
        EXAMPLES[name] = dict(n=n,k=k,minimum=cost,core_size=len(core),
                              capped=stats['capped'],rejected=stats['rejected'])
    return result


def exhaustive_suite(name,n,rankmask,oracle,weight_values,statuses,ks):
    for k in ks:
        for v in statuses:
            feas = feasible_masks(n,rankmask,k,v)
            for j,weights in enumerate(product(weight_values,repeat=n)):
                run_case(name,weights,oracle,rankmask,k,v,feas,
                         iterables=(j==0),core_property=(j % 127 == 0))


def main():
    start = time.time()
    # n=9 > rank-2,k-1 threshold 7. Exhaust all binary weights.
    for label,labels in [
        ('uniform_rank2_n9',list(range(9))),
        ('cheap_parallel_rank2_n9',[0,0,1,2,3,4,5,6,7]),
        ('loops_parallel_rank2_n9',[None,0,0,1,1,2,3,4,5])]:
        rm,ro = make_rank('classes',9,labels=labels,trunc=2)
        exhaustive_suite(label,9,rm,ro,[0,1],[(1<<9)-1,0,0b101011011],[0,1])
    # All vulnerability assignments on one fixed set of nonuniform weights.
    rm,ro = make_rank('classes',9,labels=[0,0,1,2,3,4,5,6,7],trunc=2)
    for v in range(1<<9):
        for k in (0,1,2):
            run_case('all_statuses_rank2_n9',(0,0,1,1,2,3,4,5,6),ro,rm,k,v,
                     core_property=(v % 31 == 0))
    # Genuine nonrepresentable matroid example: the Vamos matroid V8.
    rm,ro = make_rank('sparse_paving',8,r=4,forbidden=[
        (0,1,2,3),(0,1,4,5),(0,1,6,7),(2,3,4,5),(2,3,6,7)])
    exhaustive_suite('vamos_all_binary_weights',8,rm,ro,[0,1],[(1<<8)-1,0,0b10110101],[0,1])
    # Rank-3 sparse paving and a graphic matroid, independently specified oracles.
    rm,ro = make_rank('sparse_paving',9,r=3,forbidden=[
        (0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,7),(2,3,8),(2,4,6)])
    exhaustive_suite('sparse_paving_rank3',9,rm,ro,[0,1],[(1<<9)-1,0,0b101011011],[0,1])
    rm,ro = make_rank('graphic',9,vertices=4,edges=[
        (0,1),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(1,1),(2,3)])
    exhaustive_suite('graphic_rank3',9,rm,ro,[0,1],[(1<<9)-1,0,0b101011011],[0,1])
    # Rank-zero, rank-one, no elements, all loops, k much larger than n.
    for n,labels in [(0,[]),(3,[None]*3),(5,[0]*5),(6,[None,0,0,0,0,None])]:
        rm,ro=make_rank('classes',n,labels=labels,trunc=1)
        for k in (0,1,2,9):
            for v in (0,(1<<n)-1,((1<<n)-1)//3):
                for weights in [tuple([0]*n),tuple(range(n)),tuple([2]*n)]:
                    run_case('edge_cases',weights,ro,rm,k,v,iterables=True,core_property=True)
    # Explicit threshold+rejection example: only one cheap parallel representative
    # is kept by the uniform scan, while its backup must return through recursion.
    labels=[0,0,1,1,2,3,4,5,6,7]
    weights=(0,0,1,1,9,10,11,12,13,14)
    rm,ro=make_rank('classes',len(labels),labels=labels,trunc=2)
    result=run_case('adversarial_rejected_cheap_cap',weights,ro,rm,1,(1<<len(labels))-1,
                    core_property=True,iterables=True)
    assert result[3]['capped'] and result[3]['rejected']
    assert result[0]==2 and {0,1,2,3} <= result[2]
    # Returning only the capped uniform scan would cost 10, not 2.
    assert result[1] == frozenset({0,1,2,3})
    # Exact n=9 example requested alongside the text: two cheap pairs and
    # five singleton directions. The old uniform-only result costs 101.
    labels=[0,0,1,1,2,3,4,5,6]
    weights=(0,0,1,1,100,100,100,100,100)
    rm,ro=make_rank('classes',9,labels=labels,trunc=2)
    result=run_case('adversarial_n9_old_stop_101',weights,ro,rm,1,(1<<9)-1,
                    core_property=True,iterables=True)
    assert result[3]['capped'] and result[3]['rejected']
    assert result[0]==2 and result[1]==frozenset({0,1,2,3})
    old_uniform=frozenset({0,2,4,5,6,7,8})
    old_best=min(sum(weights[i] for i in members(9,b))
                 for b in feasible_masks(9,rm,1,(1<<9)-1)
                 if set(members(9,b)) <= old_uniform)
    assert old_best==101
    EXAMPLES['adversarial_n9_old_stop_101']['old_uniform_only_minimum']=old_best
    # Higher-rank threshold+rejection, exact over n=10 and k=0.
    labels=[0,0,1,2,3,4,5,6,7,8]
    weights=(0,0,1,2,3,4,5,6,7,8)
    rm,ro=make_rank('classes',10,labels=labels,trunc=3)
    result=run_case('adversarial_rank3_cap',weights,ro,rm,0,(1<<10)-1,
                    core_property=True,iterables=True)
    assert result[3]['capped'] and result[3]['rejected']
    # Every rank-2 structure below has exact checking over all subsets, with
    # wider weights including giant integers and randomized zero-heavy cases.
    rng=random.Random(20261104)
    for t in range(200):
        n=rng.randrange(1,11)
        labels=[rng.choice([None]+list(range(max(1,n//2)))) for _ in range(n)]
        trunc=rng.randrange(1,4)
        rm,ro=make_rank('classes',n,labels=labels,trunc=trunc)
        weights=tuple(rng.choice([0,0,1,2,11,10**50]) for _ in range(n))
        run_case('random_truncated_partition',weights,ro,rm,rng.randrange(4),rng.randrange(1<<n),
                 iterables=(t%11==0),core_property=True)
    # Public validation behavior and direct weighted_core generator call.
    rm,ro=make_rank('classes',3,labels=[0,1,2],trunc=2)
    for bad_k in (-1,0.5,True):
        try: solve([0,1,2],ro,bad_k)
        except ValueError: pass
        else: raise AssertionError(('bad k accepted',bad_k))
    for bad_w in ([0,-1,2],[0,1.5,2],[0,True,2]):
        try: solve(bad_w,ro,0)
        except ValueError: pass
        else: raise AssertionError(('bad weight accepted',bad_w))
    try: solve([0,1,2],ro,0,[9])
    except ValueError: pass
    else: raise AssertionError('bad vulnerable element accepted')
    c1,_=weighted_core((x for x in [0,1,2]),ro,0,(x for x in [0,2]))
    c2,_=weighted_core([0,1,2],ro,0,[0,2])
    assert c1==c2
    report=dict(counts=COUNTS,suites=SUITES,adversarial_examples=EXAMPLES,
                elapsed_seconds=round(time.time()-start,3),
                validation_cases=8,
                note='Finite verification is not proof. Weight/status exhaustiveness is per explicitly listed suite, not over all matroids.')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
