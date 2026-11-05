"""Exact subcubic realization of cotransversal robustness and disjoint quotas."""
from collections import Counter
from dataclasses import dataclass
from bounded_core_solver import source_components
from hardness_gadgets import protect_coloops
from row_splitting import dual_rank_oracle, full_rank_presentation
from scc_cover import pad_resource_network


@dataclass
class QuotaInstance:
    rows: tuple
    presentation: tuple
    elements: int
    groups: tuple
    representatives: tuple
    source_cycles: tuple
    split_rows: tuple
    split_size: int
    steps: tuple
    projection: tuple
    core_vertices: int
    core_edges: int
    leaves: int


def port_presentation(presentation, n):
    """Full-rank rows with binary root ports and a retained full matching.

    Each old row becomes a tree of rows. Its root has size at most two,
    and every other row has size three. A new column joins adjacent rows.
    The row-splitting reversal log always removes the last row and column.
    """
    original = full_rank_presentation(presentation, n)
    full = (1 << n) - 1
    covered = 0
    for row in original:
        covered |= row
    if covered != full:
        raise ValueError('The dual transversal matroid must be loopless.')
    rows, matched = list(original), [-1] * len(original)
    right_mate = [-1] * n
    for start in range(len(rows)):
        queue, seen, predecessor, free = [start], {start}, {}, None
        for left in queue:
            mask = rows[left]
            while mask:
                bit = mask & -mask
                mask ^= bit
                right = bit.bit_length() - 1
                if right in predecessor:
                    continue
                predecessor[right] = left
                other = right_mate[right]
                if other < 0:
                    free = right
                    break
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
            if free is not None:
                break
        assert free is not None
        right = free
        while right >= 0:
            left = predecessor[right]
            old = matched[left]
            matched[left], right_mate[right] = right, left
            right = old
    steps, size = [], n

    def split(i, part):
        nonlocal size
        old = rows[i]
        assert part & old == part and part and part != old
        steps.append((i, old))
        rows[i] = part | (1 << size)
        rows.append((old ^ part) | (1 << size))
        matched.append(size)
        size += 1
        return len(rows) - 1

    for start in range(len(original)):
        if rows[start].bit_count() <= 2:
            continue
        i = split(start, 1 << matched[start])
        while rows[i].bit_count() > 3:
            rest = rows[i] ^ (1 << matched[i])
            i = split(i, (1 << matched[i]) | (rest & -rest))
    return tuple(original), tuple(rows), size, tuple(steps), tuple(matched)


def quota_cotransversal_graph(presentation, n, groups, length=1):
    """Preserve robust spanning plus quota two in disjoint original subsets.

    Groups need not cover the ground set, but at least one is required.
    The target is weakly connected, loopless, bipartite, and has all
    undirected degrees in {2,3}; all locations are individually avoidable.
    Odd length multiplies the core girth and sets completion-cycle length.
    """
    groups = tuple(tuple(group) for group in groups)
    grouped = [e for group in groups for e in group]
    if (n < 2 or not groups or any(len(group) < 2 for group in groups) or
            len(grouped) != len(set(grouped)) or
            any(e < 0 or e >= n for e in grouped)):
        raise ValueError('Use nonempty disjoint quota-two groups on n elements.')
    original, split_rows, size, steps, matched = port_presentation(presentation, n)
    k = len(original)
    rows = [0] * (size + len(split_rows))
    projection = list(range(size)) + [-1] * len(split_rows)
    representatives = list(range(n))
    for i, row in enumerate(split_rows):
        center = size + i
        rows[center] = row
        rows[matched[i]] = 1 << center

    def state(origin=-1):
        v = len(rows)
        rows.append(0)
        projection.append(origin)
        return v

    cycles = []
    for group in groups:
        cycle = []
        for e in group:
            resource, filler = state(e), state()
            representatives[e] = resource
            rows[resource] = 1 << filler
            cycle.extend((resource, filler))
        for i, e in enumerate(group):
            predecessor = cycle[(2 * i - 1) % len(cycle)]
            rows[predecessor] = (1 << cycle[2 * i]) | (1 << e)
        cycles.append(tuple(cycle))

    def broadcast(parent, centers):
        if len(centers) == 1:
            rows[parent] |= 1 << centers[0]
            return
        x, y = state(), state()
        rows[parent] |= 1 << x
        rows[x], rows[y] = 1 << y, 1 << x
        middle = len(centers) // 2
        broadcast(y, centers[:middle])
        broadcast(y, centers[middle:])

    broadcast(representatives[groups[0][0]], list(range(size, size + k)))
    for v in range(n):
        protected = rows[v]  # The matched center stays bidirectionally adjacent.
        while True:
            incoming = [u for u, row in enumerate(rows) if row >> v & 1]
            if len(incoming) <= 3:
                break
            moved = [u for u in incoming if not protected >> u & 1][:2]
            assert len(moved) == 2
            x, y = state(v), state()
            rows[x], rows[y] = 1 << y, (1 << x) | (1 << v)
            for u in moved:
                rows[u] ^= (1 << v) | (1 << x)

    leaves = [v for v in range(n) if sum(
        bool(rows[v] >> u & 1 or rows[u] >> v & 1) for u in range(len(rows))) == 1]
    core_vertices = len(rows)
    core_edges = sum(bool(rows[u] >> v & 1 or rows[v] >> u & 1)
                     for u in range(len(rows)) for v in range(u + 1, len(rows)))
    rows, cycles, _, projection, _ = pad_resource_network(
        tuple(rows), tuple(cycles), tuple(int(e >= 0) for e in projection),
        tuple(projection), {}, length)
    rows, projection = list(rows), list(projection)
    # Balanced cycle targets are parallel to a leaf target; centers are loops.
    for v in leaves:
        cycle = [state(v if i % 2 == 0 else -1) for i in range(4 * length)]
        for i, u in enumerate(cycle):
            w = cycle[(i + 1) % len(cycle)]
            rows[u] |= 1 << w
            rows[w] |= 1 << u
        center = cycle[1]
        rows[v] |= 1 << center
        rows[center] |= 1 << v
    return QuotaInstance(tuple(rows), original, n, groups, tuple(representatives),
                         cycles, split_rows, size, steps, tuple(projection),
                         core_vertices, core_edges, len(leaves))


def quota_from_placement(instance, selected):
    """Retain all selected original representatives while reversing the map.

    Feasibility of the target placement is a precondition. No selected
    original representative is removed by source normalization, coloop
    protection or row-splitting exchanges.
    """
    selected = set(selected)
    if not selected <= set(range(len(instance.rows))):
        raise ValueError('Invalid state index.')
    initial_size = len(selected)
    for group, cycle in zip(instance.groups, instance.source_cycles):
        group = {instance.representatives[e] for e in group}
        cycle = set(cycle)
        if len(selected & cycle) < 2:
            raise ValueError('A source quota is not satisfied.')
        assert all(instance.projection[v] == -1 for v in cycle - group)
        selected.difference_update(cycle - group)
        for v in sorted(group - selected):
            if len(selected & group) >= 2:
                break
            selected.add(v)
    counts = Counter(instance.projection[v] for v in selected if instance.projection[v] >= 0)
    support = set(protect_coloops(dual_rank_oracle(instance.split_rows, instance.split_size),
                                 range(instance.split_size), counts))
    assert len(support) <= len(selected) <= initial_size
    rows, size = list(instance.split_rows), instance.split_size
    d = instance.elements - len(instance.presentation)
    for i, old in reversed(instance.steps):
        rows = rows[:-1]
        rows[i] = old
        size -= 1
        if size not in support:
            continue
        support.remove(size)
        rank = dual_rank_oracle(rows, size)

        def robust(chosen):
            return rank(chosen) == d and all(rank(chosen - {e}) == d for e in chosen)

        if not robust(support):
            replacement = next((w for w in range(size) if w not in support and
                                robust(support | {w})), None)
            if replacement is None:
                raise AssertionError('The row-splitting exchange failed.')
            support.add(replacement)
    assert len(support) <= initial_size
    assert all(len(support & set(group)) >= 2 for group in instance.groups)
    return frozenset(support)


def sparse_control_graph(rows, length=1):
    """Exact sparse realization of any feasible general directed control graph."""
    rows = tuple(rows)
    n = len(rows)
    groups = tuple(tuple(v for v in range(n) if block >> v & 1)
                   for block in source_components(rows))
    return quota_cotransversal_graph(rows, n, groups, length)
