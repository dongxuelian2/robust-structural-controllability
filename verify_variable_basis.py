"""Exact audits of OPT=n+pc(H) and polynomial solution extraction."""
from itertools import combinations
from random import Random
from variable_basis import (realize, sparse_realize, improve_cover,
                            cover_from_placement, cover_from_parallel,
                            placement_from_cover)
from hardness_gadgets import subdivide_ports
from audit_hardness import augmented_matchable


def path_partition(n, edges):
    neighbors = [0] * n
    for u, v in edges:
        neighbors[u] |= 1 << v
        neighbors[v] |= 1 << u
    paths = [None] * (1 << n)
    endpoint = [dict() for _ in paths]
    for v in range(n):
        endpoint[1 << v][v] = (v,)
    for mask in range(1, 1 << n):
        for last, path in tuple(endpoint[mask].items()):
            paths[mask] = path
            available = neighbors[last] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                v = bit.bit_length() - 1
                endpoint[mask | bit][v] = path + (v,)
    optimum = [None] * (1 << n)
    optimum[0] = ()
    for mask in range(1, 1 << n):
        first = mask & -mask
        part = mask
        while part:
            if part & first and paths[part] is not None:
                candidate = optimum[mask ^ part] + ((False, paths[part]),)
                if optimum[mask] is None or len(candidate) < len(optimum[mask]):
                    optimum[mask] = candidate
            part = (part - 1) & mask
    return optimum[-1]


def assignment_ranks(n, edges):
    """Direct assignment matroid ranks, no state-dual or graph-cycle formula."""
    m = len(edges)
    job_rows = [1 << v for v in range(n)]
    for j, (u, v) in enumerate(edges):
        job_rows[u] |= 1 << (n + j)
        job_rows[v] |= 1 << (n + j)
    reachable = {0}
    for row in job_rows:
        previous = tuple(reachable)
        for mask in previous:
            remaining = row & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                reachable.add(mask | bit)
    ranks = [mask.bit_count() if mask in reachable else 0
             for mask in range(1 << (n + m))]
    for j in range(n + m):
        for mask in range(1 << (n + m)):
            if mask >> j & 1:
                ranks[mask] = max(ranks[mask], ranks[mask ^ (1 << j)])
    return ranks


def check_cover(n, edges, cover):
    edge_set = set(edges)
    vertices = []
    for cycle, part in cover:
        vertices.extend(part)
        assert all(tuple(sorted(e)) in edge_set for e in zip(part, part[1:]))
        if cycle:
            assert len(part) >= 3 and tuple(sorted((part[-1], part[0]))) in edge_set
    assert sorted(vertices) == list(range(n))


def check(n, edges, direct=False):
    rows, edges = realize(n, edges)
    ranks = assignment_ranks(n, edges)
    cover = path_partition(n, edges)
    expected = n + len(cover)
    best, repeated_best = n + len(edges), 2 * n
    intersection = (1 << (n + len(edges))) - 1
    for mask, rank in enumerate(ranks):
        selected = {j for j in range(n + len(edges)) if mask >> j & 1}
        critical = [v for v in selected if ranks[mask ^ (1 << v)] < n] if rank == n else []
        valid = rank == n and not critical
        if rank == n:
            cost = len(selected) + len(critical)
            repeated_best = min(repeated_best, cost)
            chosen_edges = [e for j, e in enumerate(edges) if n + j in selected]
            decoded = improve_cover(n, chosen_edges)
            check_cover(n, edges, decoded)
            assert len(decoded) <= cost - n, (n, edges, mask, critical, decoded)
        if direct:
            actual = augmented_matchable(rows, selected) and all(
                augmented_matchable(rows, selected - {v}) for v in selected)
            assert valid == actual, (n, edges, selected, 'state-model equivalence')
        if valid:
            best = min(best, len(selected))
            intersection &= mask
            decoded = cover_from_placement(n, edges, selected)
            check_cover(n, edges, decoded)
            assert len(decoded) <= len(selected) - n
    assert best == repeated_best == expected, (n, edges, best, repeated_best, expected)
    if all(sum(v in e for e in edges) >= 2 for v in range(n)):
        assert intersection == 0, ('a mandatory resource', n, edges, intersection)
    lifted = placement_from_cover(n, edges, cover)
    lifted_mask = sum(1 << v for v in lifted)
    assert len(lifted) == expected
    assert ranks[lifted_mask] == n and all(ranks[lifted_mask ^ (1 << v)] == n
                                         for v in lifted)


def audit_sparse_and_padded(n, edges, rng):
    base_rows, edges = realize(n, edges)
    sparse = sparse_realize(n, edges)
    padded = subdivide_ports(sparse)
    ranks = assignment_ranks(n, edges)
    for instance in (sparse, padded):
        assert all(row.bit_count() <= 3 for row in instance.rows)
        full = (1 << len(instance.rows)) - 1
        # Use the established BFS matcher here, separately from the job oracle.
        from bounded_core_solver import matching_number
        nu = matching_number(instance.rows, full)
        assert len(instance.rows) - nu == n
        for _ in range(100):
            selected = {v for v in range(len(instance.rows)) if rng.randrange(2)}
            counts = {}
            for v in selected:
                origin = instance.projection[v]
                if origin >= 0:
                    counts[origin] = counts.get(origin, 0) + 1
            support = sum(1 << v for v in counts)
            selected_mask = sum(1 << v for v in selected)
            actual_rank = len(selected) + matching_number(
                instance.rows, full ^ selected_mask) - nu
            assert actual_rank == ranks[support], (n, edges, 'sparse rank identity')
            if ranks[support] == n and all(ranks[support ^ (1 << v)] == n
                                           for v, amount in counts.items() if amount == 1):
                decoded = cover_from_parallel(n, edges, instance, selected)
                check_cover(n, edges, decoded)
                assert n + len(decoded) <= len(selected)


if __name__ == '__main__':
    count = 0
    for n in range(2, 6):
        possible = list(combinations(range(n), 2))
        for bits in range(1, 1 << len(possible)):
            edges = tuple(e for j, e in enumerate(possible) if bits >> j & 1)
            if {v for e in edges for v in e} != set(range(n)):
                continue
            check(n, edges, direct=(n <= 3))
            count += 1
        print('All isolate-free source graphs through n=%d: %d PASS' % (n, count),
              flush=True)
    rng = Random(20261104)
    random_count = 0
    for _ in range(50):
        n = 6
        edges = tuple(e for e in combinations(range(n), 2) if rng.randrange(3) == 0)
        if {v for e in edges for v in e} == set(range(n)):
            check(n, edges)
            random_count += 1
    print('Additional n=6 fixed-seed sparse graphs:', random_count, 'PASS', flush=True)
    for n, edges in ((3, ((0, 1), (0, 2), (1, 2))),
                     (4, tuple(combinations(range(4), 2))),
                     (5, tuple((0, v) for v in range(1, 5)))):
        audit_sparse_and_padded(n, edges, rng)
    print('Sparse/padded realization: 600 full-state rank queries and reverse maps PASS')
