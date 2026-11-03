"""Cubic-time factor-two approximation for single dedicated-actuator failure.

Rows are outgoing adjacency bitmasks. Returns placement and verifiable statistics.
No repeated actuators; infeasible graphs return None.
"""
from collections import deque
from bounded_core_solver import source_components


def solve(rows):
    rows = tuple(rows)
    n = len(rows)
    if n == 0 or any(not isinstance(row, int) or row < 0 or row >= 1 << n for row in rows):
        raise ValueError('Expected a nonempty sequence of nonnegative n-bit rows.')
    sources = source_components(rows)
    if any(c.bit_count() < 2 for c in sources):
        return None
    # s=0, left states=1..n, right states=n+1..2n,
    # group nodes=2n+1..2n+ell, then sink.
    sink = 2*n + len(sources) + 1
    cap = [[0]*(sink+1) for _ in range(sink+1)]
    group = {}
    for j, component in enumerate(sources):
        for i in range(n):
            if component >> i & 1: group[i] = 2*n+1+j
        cap[2*n+1+j][sink] = component.bit_count()-2
    for i,row in enumerate(rows):
        cap[0][1+i] = 1
        cap[n+1+i][group.get(i,sink)] = 1
        for j in range(n):
            if row >> j & 1: cap[1+i][n+1+j] = 1
    residual = [row[:] for row in cap]
    flow = 0
    while True:
        parent = [-1]*(sink+1); parent[0]=0; queue=deque([0])
        while queue and parent[sink] < 0:
            u = queue.popleft()
            for w,c in enumerate(residual[u]):
                if c and parent[w] < 0:
                    parent[w]=u; queue.append(w)
        if parent[sink] < 0: break
        # Every source-to-sink augmenting path has residual bottleneck one.
        w = sink
        while w:
            u=parent[w]; residual[u][w]-=1; residual[w][u]+=1; w=u
        flow += 1
    left_mate = [-1]*n
    covered = 0
    for i in range(n):
        for j in range(n):
            if cap[1+i][n+1+j] and residual[1+i][n+1+j]==0:
                left_mate[i] = j
                covered |= 1<<j
    assert covered.bit_count() == flow
    base = ((1<<n)-1) ^ covered
    incoming = [[i for i,row in enumerate(rows) if row>>j&1] for j in range(n)]
    critical=[]; repairs=[]
    for b in range(n):
        if not base>>b&1: continue
        # Is there an augmenting path that adds right vertex b to this matching?
        seen_right={b}; seen_left=set(); queue=deque([b]); augmentable=False
        while queue and not augmentable:
            right=queue.popleft()
            for left in incoming[right]:
                if left in seen_left: continue
                seen_left.add(left)
                matched=left_mate[left]
                if matched < 0:
                    augmentable=True; break
                if matched not in seen_right:
                    seen_right.add(matched); queue.append(matched)
        if not augmentable:
            critical.append(b)
            assert incoming[b], 'Feasibility guarantees every state has an incoming edge.'
            replacement=left_mate[incoming[b][0]]
            assert replacement >= 0 and covered>>replacement&1
            repairs.append(replacement)
    chosen=base
    for e in repairs: chosen |= 1<<e
    return {
        'actuators': tuple(i for i in range(n) if chosen>>i&1),
        'relaxation_actuators': tuple(i for i in range(n) if base>>i&1),
        'relaxation_value': base.bit_count(),
        'critical_states': tuple(critical),
        'repair_states': tuple(sorted(set(repairs))),
        'matching': tuple((i,j) for i,j in enumerate(left_mate) if j>=0),
    }
