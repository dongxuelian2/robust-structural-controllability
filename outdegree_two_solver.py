"""Exact unit-cost 1-RMSCP-d for state outdegree at most two; Python 3.10+."""
from bounded_core_solver import source_components


def tree_blocks(rows):
    """Tree components of the multigraph whose edges are nonempty row supports.

    A one-neighbor row becomes a loop; distinct rows may become parallel edges.
    Return their state-vertex sets as bitmasks.
    """
    n = len(rows)
    parent = list(range(n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    supports = []
    for row in rows:
        support = [j for j in range(n) if row >> j & 1]
        if len(support) > 2:
            raise ValueError('The state outdegree must be at most two.')
        if support:
            supports.append(support)
            if len(support) == 2:
                a, b = map(find, support)
                parent[a] = b
    blocks, edge_count = {}, {}
    for v in range(n):
        r = find(v)
        blocks[r] = blocks.get(r, 0) | 1 << v
    for support in supports:
        r = find(support[0])
        edge_count[r] = edge_count.get(r, 0) + 1
    return tuple(mask for r, mask in blocks.items()
                 if edge_count.get(r, 0) == mask.bit_count() - 1)


def solve(rows):
    """Return an optimum tuple of distinct actuator states, or None.

    Matching robustness supplies one partition of two-element quotas; source
    SCCs supply a second partition. A capacitated bipartite matching maximizes
    how many of the two sets of quotas can be met by shared actuator states.
    """
    rows = tuple(rows)
    n = len(rows)
    full = (1 << n) - 1
    if not n or any(not isinstance(r, int) or r < 0 or r & ~full for r in rows):
        raise ValueError('Use a nonempty sequence of nonnegative n-bit rows.')
    trees = tree_blocks(rows)
    sources = source_components(rows)
    if any(c.bit_count() < 2 for c in sources):
        return None
    if any(c.bit_count() < 2 for c in trees):
        raise AssertionError('An isolated right state has a singleton source SCC.')

    # A state node gives each intersection edge capacity one even when both
    # of its quota-group endpoints have capacity two.
    a_count, b_count = len(trees), len(sources)
    start = 0
    tree_offset = 1
    state_offset = tree_offset + a_count
    source_offset = state_offset + n
    sink = source_offset + b_count
    adjacency = [[] for _ in range(sink + 1)]
    capacities = {}

    def add(a, b, capacity):
        adjacency[a].append(b)
        adjacency[b].append(a)
        capacities[a, b] = capacity
        capacities[b, a] = 0

    for i, block in enumerate(trees):
        add(start, tree_offset + i, 2)
        for v in range(n):
            if block >> v & 1:
                add(tree_offset + i, state_offset + v, 1)
    for j, block in enumerate(sources):
        add(source_offset + j, sink, 2)
        for v in range(n):
            if block >> v & 1:
                add(state_offset + v, source_offset + j, 1)
    overlap = 0
    while True:
        predecessor = {start: None}
        queue = [start]
        for u in queue:
            for v in adjacency[u]:
                if v not in predecessor and capacities[u, v]:
                    predecessor[v] = u
                    queue.append(v)
            if sink in predecessor:
                break
        if sink not in predecessor:
            break
        v = sink
        while v != start:
            u = predecessor[v]
            capacities[u, v] -= 1
            capacities[v, u] += 1
            v = u
        overlap += 1

    selected = 0
    for i, block in enumerate(trees):
        for v in range(n):
            if block >> v & 1 and capacities.get((state_offset + v,
                                                tree_offset + i), 0):
                selected |= 1 << v
    for block in (*trees, *sources):
        while (selected & block).bit_count() < 2:
            available = block & ~selected
            selected |= available & -available
    assert selected.bit_count() == 2 * (a_count + b_count) - overlap
    return tuple(v for v in range(n) if selected >> v & 1)
