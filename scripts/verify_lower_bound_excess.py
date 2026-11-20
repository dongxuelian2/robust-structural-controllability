"""Exact checks for the Hamilton-path/absolute-lower-bound reduction.

Python 3.10+, standard library only. Actual state/input matching and reachability
are checked independently of the matroid proof. Exhaustive checks are finite
supplements, not the NP-hardness proof.
"""
from itertools import combinations, permutations
from random import Random
from collections import Counter


def matching_size(rows, nright):
    """Independent augmenting-path matching, rows are integer neighborhoods."""
    mate = [-1] * nright
    def add(row, seen):
        mask = rows[row]
        while mask:
            bit = mask & -mask
            mask -= bit
            v = bit.bit_length() - 1
            if seen[v]:
                continue
            seen[v] = True
            if mate[v] < 0 or add(mate[v], seen):
                mate[v] = row
                return True
        return False
    return sum(add(i, [False] * nright) for i in range(len(rows)))


def assignment_robust(resources, n):
    if matching_size(resources, n) != n:
        return False
    return all(matching_size(resources[:i] + resources[i+1:], n) == n
               for i in range(len(resources)))


def hamilton_path(n, edges, s=0, t=None):
    t = n-1 if t is None else t
    es = {frozenset(e) for e in edges}
    for middle in permutations(v for v in range(n) if v not in (s, t)):
        path = (s,) + middle + (t,)
        if all(frozenset((u,v)) in es for u,v in zip(path,path[1:])):
            return path
    return None


def build_control(n, edges, s=0, t=None):
    t = n-1 if t is None else t
    resources = [(1 << u) | (1 << v) for u,v in edges] + [1 << s, 1 << t]
    m = len(resources)
    # resources [0,m), centers [m,2m), jobs [2m,2m+n), source A,B,FA,FB.
    A, B, FA, FB = range(2*m+n, 2*m+n+4)
    rows = [0] * (2*m+n+4)
    for e, nb in enumerate(resources):
        rows[e] = 1 << (m+e)
        rows[m+e] = (1 << e) | (nb << (2*m))
    rows[A] = (1 << FA) | (((1 << m)-1) << m)
    rows[FA] = (1 << B) | (1 << (m-1))
    rows[B] = 1 << FB
    rows[FB] = (1 << A) | (1 << (m-2))
    return tuple(rows), resources, (A,B,FA,FB)


def actual_control_robust(rows, selected):
    """Use full state-plus-input graph and direct accessibility, no rank formula."""
    N = len(rows)
    full = (1 << N)-1
    for failed in (-1,) + tuple(selected):
        survivors = [v for v in selected if v != failed]
        reach = sum(1 << v for v in survivors)
        old = -1
        while old != reach:
            old = reach
            for v in range(N):
                if reach >> v & 1:
                    reach |= rows[v]
        if reach != full:
            return False
        if matching_size(rows + tuple(1 << v for v in survivors), N) != N:
            return False
    return True


def source_sccs(rows):
    N = len(rows)
    reach = []
    for u in range(N):
        r = 1 << u
        old = -1
        while r != old:
            old = r
            for v in range(N):
                if r >> v & 1:
                    r |= rows[v]
        reach.append(r)
    remain = set(range(N))
    comps = []
    while remain:
        u = min(remain)
        C = {v for v in remain if reach[u] >> v & 1 and reach[v] >> u & 1}
        remain -= C
        comps.append(C)
    return [C for C in comps if not any(rows[u] >> v & 1
            for u in range(N) if u not in C for v in C)]


def project_resources(selected, resources, source, n):
    m = len(resources)
    A,B,FA,FB = source
    out = []
    for v in selected:
        if v < m:
            out.append(resources[v])
        elif 2*m <= v < 2*m+n:
            out.append(1 << (v-2*m))
        elif v == A:
            out.append(resources[-2])
        elif v == B:
            out.append(resources[-1])
        else:
            out.append(0)
    return tuple(out)


def clone_vertex(n, edges, v=0):
    old = [u for u in range(n) if u != v]
    ix = {u:i+1 for i,u in enumerate(old)}
    s,t = 0,n
    out = []
    for a,b in edges:
        if a == v or b == v:
            other = b if a == v else a
            out.extend(((s,ix[other]), (t,ix[other])))
        else:
            out.append((ix[a],ix[b]))
    return n+1, sorted(out)


def main():
    counts = Counter()
    # Exhaust every graph through five vertices and every anchored threshold
    # support, including every possible completion singleton.
    for n in range(2,6):
        pairs = list(combinations(range(n),2))
        for mask in range(1 << len(pairs)):
            edges = [e for i,e in enumerate(pairs) if mask >> i & 1]
            options = tuple((1 << u)|(1 << v) for u,v in edges) + tuple(1 << u for u in range(n))
            found = False
            for ids in combinations(range(len(options)),n-1):
                selected = (1, 1 << (n-1)) + tuple(options[i] for i in ids)
                good = assignment_robust(selected,n)
                counts['anchored_assignment_supports'] += 1
                if good:
                    found = True
                    assert all(i < len(edges) for i in ids), (n,edges,ids)
                    chosen = [edges[i] for i in ids]
                    assert hamilton_path(n,chosen), (n,edges,ids)
            assert found == bool(hamilton_path(n,edges)), (n,edges)
            counts['all_small_H_graphs'] += 1
    print('Source threshold checks:', dict(counts), flush=True)

    # Exhaust actual control placements on every non-isolated H through order
    # three and on ten fixed-seed order-four graphs. No matroid-loop filtering.
    rng = Random(20261106)
    fixtures = []
    for n in (2,3):
        pairs = list(combinations(range(n),2))
        for mask in range(1 << len(pairs)):
            edges = [e for i,e in enumerate(pairs) if mask >> i & 1]
            if all(any(v in e for e in edges) for v in range(n)):
                fixtures.append((n,edges))
    pairs = list(combinations(range(4),2))
    masks = list(range(1 << len(pairs)))
    rng.shuffle(masks)
    for mask in masks:
        edges = [e for i,e in enumerate(pairs) if mask >> i & 1]
        if all(any(v in e for e in edges) for v in range(4)):
            fixtures.append((4,edges))
            if len(fixtures) == 15:
                break
    for n,edges in fixtures:
        rows,resources,source = build_control(n,edges)
        N = len(rows)
        assert N - matching_size(rows,N) == n
        assert source_sccs(rows) == [set(source)]
        found = False
        for selected in combinations(range(N),n+1):
            if len(set(selected) & set(source)) < 2:
                continue
            good = actual_control_robust(rows,selected)
            expected = (len(set(selected) & set(source)) >= 2 and
                        assignment_robust(project_resources(selected,resources,source,n),n))
            assert good == expected, (n,edges,selected)
            counts['actual_control_threshold_placements'] += 1
            if good:
                found = True
                assert source[0] in selected and source[1] in selected
        assert found == bool(hamilton_path(n,edges)), (n,edges)
        counts['exhaustive_control_fixtures'] += 1
    print('Control threshold checks:', dict(counts), flush=True)

    # Hamiltonian-cycle reductions, including positive and negative fixtures.
    cube = [(u,u^(1<<b)) for u in range(8) for b in range(3) if u < (u^(1<<b))]
    petersen = [(i,(i+1)%5) for i in range(5)] + [(i,i+5) for i in range(5)] + [(i+5,(i+2)%5+5) for i in range(5)]
    examples = [
        (4,list(combinations(range(4),2)),True),
        (6,[(u,v) for u in range(3) for v in range(3,6)],True),
        (8,cube,True),
        (10,petersen,False),
        (5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],False),
    ]
    for old_n,old_edges,yes in examples:
        n,edges = clone_vertex(old_n,old_edges)
        path = hamilton_path(n,edges)
        assert bool(path) == yes
        rows,resources,source = build_control(n,edges)
        N = len(rows)
        assert N-matching_size(rows,N) == n
        assert source_sccs(rows) == [set(source)]
        assert all(not (rows[u] >> u & 1) for u in range(N))
        m = len(resources)
        X = set(range(m)) | set(range(2*m,2*m+n)) | set(source[:2])
        assert all((u in X) != (v in X) for u in range(N) for v in range(N) if rows[u] >> v & 1)
        if path:
            index = {frozenset(e):i for i,e in enumerate(edges)}
            selected = tuple(index[frozenset((a,b))] for a,b in zip(path,path[1:])) + source[:2]
            assert len(selected) == n+1
            assert actual_control_robust(rows,selected)
            counts['positive_lift_fixtures'] += 1
        # All pairs of omitted state columns match; verify directly, not by
        # relying on their displayed neighborhoods. This plus source size four
        # certifies every individual location is avoidable under every failure.
        for u,v in combinations(range(N),2):
            restricted = tuple(row & ((1<<u)|(1<<v)) for row in rows)
            assert matching_size(restricted,N) == 2
            counts['omitted_state_pair_checks'] += 1
        for _ in range(80):
            selected = tuple(rng.sample(range(N),rng.randrange(min(N,n+5)+1)))
            good = actual_control_robust(rows,selected)
            expected = (len(set(selected)&set(source)) >= 2 and
                        assignment_robust(project_resources(selected,resources,source,n),n))
            assert good == expected
            counts['arbitrary_control_placements'] += 1
        counts['cloned_Hamilton_fixtures'] += 1
    print('PASS', dict(sorted(counts.items())))


if __name__ == '__main__':
    main()
