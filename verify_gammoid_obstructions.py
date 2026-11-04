"""Exact certificates beyond private-basis assignment models; standard library."""
from collections import deque
from itertools import combinations, product
from random import Random
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number
from hardness_gadgets import cotransversal_graph, protect_coloops, subdivide_ports


def transversal_ranks(rows, n):
    reachable = {0}
    for row in rows:
        for mask in tuple(reachable):
            remaining = row & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                reachable.add(mask | bit)
    ranks = [mask.bit_count() if mask in reachable else 0
             for mask in range(1 << n)]
    for v in range(n):
        for mask in range(1 << n):
            if mask >> v & 1:
                ranks[mask] = max(ranks[mask], ranks[mask ^ (1 << v)])
    return ranks


def dual_ranks(rows, n):
    ranks = transversal_ranks(rows, n)
    full = (1 << n) - 1
    assert ranks[full] == len(rows)
    return [mask.bit_count() + ranks[full ^ mask] - len(rows)
            for mask in range(1 << n)]


def robust(ranks, mask):
    return ranks[mask] == ranks[-1] and all(
        ranks[mask ^ (1 << v)] == ranks[-1]
        for v in range(len(ranks).bit_length() - 1) if mask >> v & 1)


def state_rows(presentation, n):
    """Bidirected incidence realization: right centers are a free summand."""
    result = [0] * (n + len(presentation))
    for j, row in enumerate(presentation):
        result[n + j] = row
        for v in range(n):
            if row >> v & 1:
                result[v] |= 1 << (n + j)
    return tuple(result)


def cyclic_flat(ranks, mask):
    n = len(ranks).bit_length() - 1
    r = ranks[mask]
    return all(ranks[mask | (1 << v)] > r for v in range(n) if not mask >> v & 1) and all(
        ranks[mask ^ (1 << v)] == r for v in range(n) if mask >> v & 1)


def component_count(ranks):
    n = len(ranks).bit_length() - 1
    d = ranks[-1]
    basis = next(mask for mask, r in enumerate(ranks) if r == d and mask.bit_count() == d)
    adjacency = [0] * n
    for v in range(n):
        if basis >> v & 1:
            continue
        for b in range(n):
            if basis >> b & 1 and ranks[(basis ^ (1 << b)) | (1 << v)] == d:
                adjacency[v] |= 1 << b
                adjacency[b] |= 1 << v
    unseen, count = (1 << n) - 1, 0
    while unseen:
        reached = unseen & -unseen
        while True:
            extended = reached
            for v in range(n):
                if reached >> v & 1:
                    extended |= adjacency[v]
            if extended == reached:
                break
            reached = extended
        unseen &= ~reached
        count += 1
    return count


def audit_nontransversal():
    presentation = (1033, 4612, 2084, 1064, 560, 352, 448, 4114)
    n = 13
    ranks = dual_ranks(presentation, n)
    tr = transversal_ranks(presentation, n)
    assert ranks[-1] == 5
    assert all(ranks[(1 << u) | (1 << v)] == tr[(1 << u) | (1 << v)] == 2
               for u, v in combinations(range(n), 2))
    assert component_count(ranks) == 1
    a, b, c = 480, 4114, 4612
    assert all(cyclic_flat(ranks, mask) for mask in (a, b, c))
    assert [ranks[x] for x in (a, b, c, a | b, a | c, b | c, a | b | c, a & b & c)] == [
        2, 2, 2, 4, 4, 3, 4, 0]
    # Mason's necessary transversal inequality would require 0 <= -1.
    optimum = min(mask.bit_count() for mask in range(1 << n) if robust(ranks, mask))
    assert optimum == 7 and robust(ranks, 223)
    states = state_rows(presentation, n)
    for subset in combinations(range(n), 6):
        selected = set(subset)
        assert not all(augmented_matchable(states, selected - {v}) for v in selected)
    selected = {v for v in range(n) if 223 >> v & 1}
    assert all(augmented_matchable(states, selected - {v}) for v in selected)

    # Split the sole degree-four target, moving its incidences in rows 2,3.
    sparse = list(presentation)
    for j in (2, 3):
        sparse[j] ^= (1 << 5) | (1 << 13)
    sparse.append((1 << 5) | (1 << 13))
    sparse_ranks = dual_ranks(sparse, 14)
    for mask, rank in enumerate(sparse_ranks):
        support = mask & ((1 << 13) - 1)
        if mask >> 13 & 1:
            support |= 1 << 5
        assert rank == ranks[support]
    sparse_states = state_rows(sparse, 14)
    assert len(sparse_states) == 23 and max(row.bit_count() for row in sparse_states) == 3
    assert min(mask.bit_count() for mask in range(1 << 14)
               if robust(sparse_ranks, mask)) == 7
    assert all(augmented_matchable(sparse_states, selected - {v}) for v in selected)
    print('Connected simple/cosimple nontransversal prototype: rank 5, OPT 7;'
          ' 1716 direct six-input exclusions; 16384 split ranks PASS', flush=True)


def two_hub_network(n, arcs):
    arcs = tuple(arcs)
    count = 3 * n + 2 + len(arcs)
    rows = [0] * count
    for v in range(n):
        rows[n + 2 + v] = (1 << v) | (1 << n)
        rows[2 * n + 2 + v] = (1 << v) | (1 << (n + 1))
    for j, (u, v) in enumerate(arcs):
        rows[3 * n + 2 + j] = (1 << (n + 2 + u)) | (1 << (2 * n + 2 + v))
    return tuple(rows), frozenset(range(n + 2))


def linkage_rank(rows, sinks, selected):
    """Independent integer vertex-capacitated max-flow oracle."""
    n = len(rows)
    source, target = 2 * n, 2 * n + 1
    residual = [dict() for _ in range(2 * n + 2)]

    def edge(u, v):
        residual[u][v] = 1
        residual[v].setdefault(u, 0)

    for v in range(n):
        edge(2 * v, 2 * v + 1)
        if v in selected:
            edge(source, 2 * v)
        if v in sinks:
            edge(2 * v + 1, target)
        for w in range(n):
            if rows[v] >> w & 1:
                edge(2 * v + 1, 2 * w)
    flow = 0
    while True:
        parent = {source: None}
        queue = deque([source])
        while queue and target not in parent:
            u = queue.popleft()
            for v, capacity in residual[u].items():
                if capacity and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if target not in parent:
            return flow
        v = target
        while v != source:
            u = parent[v]
            residual[u][v] -= 1
            residual[v][u] += 1
            v = u
        flow += 1


def audit_directed_comb(k):
    """Strongly connected source; exact path-cover 2k-1, robust target 4k."""
    n, z = 3 * k + 1, 3 * k
    leaves = [k + 2 * v + i for v in range(k) for i in (0, 1)]
    arcs = ([(v, v + 1) for v in range(k - 1)] +
            [(v, k + 2 * v + i) for v in range(k) for i in (0, 1)] +
            [(leaf, z) for leaf in leaves] + [(z, 0)])
    for reverse in (False, True):
        reached, queue = {0}, [0]
        for u in queue:
            for a, b in arcs:
                if reverse:
                    a, b = b, a
                if a == u and b not in reached:
                    reached.add(b)
                    queue.append(b)
        assert len(reached) == n
    rows, sinks = two_hub_network(n, arcs)
    index = {arc: 3 * n + 2 + j for j, arc in enumerate(arcs)}
    selected = {n + 2 + leaf for leaf in leaves}
    selected |= {index[v, k + 2 * v + i] for v in range(k) for i in (0, 1)}
    for v in (0, 1):
        leaf = k + 2 * v
        selected.remove(n + 2 + leaf)
        selected.add(index[leaf, z])
    assert len(selected) == 4 * k
    for failed in selected:
        assert linkage_rank(rows, sinks, selected - {failed}) == n + 2

    # Every leaf sink has five possible starting resources. These groups
    # are disjoint, proving at least two actuators per leaf are necessary.
    groups = [{leaf, n + 2 + leaf, 2 * n + 2 + leaf,
               index[v, leaf], index[leaf, z]}
              for v in range(k) for leaf in (k + 2 * v, k + 2 * v + 1)]
    assert len(set.union(*groups)) == sum(map(len, groups))
    for leaf, group in zip(leaves, groups):
        reachable = {leaf}
        while True:
            extra = {v for v, row in enumerate(rows)
                     if any(row >> w & 1 for w in reachable)}
            if extra <= reachable:
                break
            reachable |= extra
        assert reachable == group
    special = [leaves[0], z] + list(range(k)) + [leaves[-2]]
    cover = [special] + [[leaf] for leaf in leaves if leaf not in special]
    assert len(cover) == 2 * k - 1
    assert sorted(v for path in cover for v in path) == list(range(n))
    assert all((u, v) in arcs for path in cover for u, v in zip(path, path[1:]))
    print('Strong directed comb k =', k, ': pc =', 2 * k - 1,
          '; target OPT =', 4 * k, '; proposed n+pc =', 5 * k, 'PASS', flush=True)


def audit_general_realization():
    rng = Random(20261104)
    count = queries = reverse = 0

    def check(presentation, n, padding=False):
        nonlocal count, queries, reverse
        instance = cotransversal_graph(presentation, n)
        r = len(presentation)
        original = dual_ranks(presentation, n)
        d = original[-1]
        assert len(instance.rows) <= 11 * r - n
        assert sum(row.bit_count() for row in instance.rows) // 2 <= 14 * r - 2 * n - 1
        if padding:
            instance = subdivide_ports(instance)
        rows, projection = instance.rows, instance.projection
        size, full = len(rows), (1 << len(rows)) - 1
        assert max(row.bit_count() for row in rows) <= 3
        color, queue = {0: 0}, [0]
        for u in queue:
            for v in range(size):
                if not rows[u] >> v & 1:
                    continue
                assert rows[v] >> u & 1
                if v not in color:
                    color[v] = 1 - color[u]
                    queue.append(v)
                assert color[v] != color[u]
        assert len(color) == size
        nu = matching_number(rows, full)
        assert size - nu == d

        def actual(selected):
            nonlocal queries
            queries += 1
            mask = sum(1 << v for v in selected)
            return len(selected) + matching_number(rows, full ^ mask) - nu

        supports = [set(v for v in range(n) if mask >> v & 1)
                    for mask in range(1 << n)] if not padding else []
        supports += [{v for v in range(size) if rng.randrange(2)} for _ in range(2)]
        for selected in supports:
            counts = {}
            for v in selected:
                if projection[v] >= 0:
                    counts[projection[v]] = counts.get(projection[v], 0) + 1
            support = sum(1 << v for v in counts)
            assert actual(selected) == original[support]
            valid = original[support] == d and all(
                original[support ^ (1 << v)] == d for v, amount in counts.items() if amount == 1)
            if valid:
                for v in selected:
                    assert actual(selected - {v}) == d
                normalized = protect_coloops(
                    lambda states: original[sum(1 << v for v in states)], range(n), counts)
                assert len(normalized) <= len(selected)
                assert robust(original, sum(1 << v for v in normalized))
                reverse += 1
        count += 1

    for n in range(2, 5):
        full = (1 << n) - 1
        choices = [mask for mask in range(1, full + 1) if mask.bit_count() <= 3]
        for r in range(1, min(n, 4)):
            for presentation in product(choices, repeat=r):
                covered = 0
                for row in presentation:
                    covered |= row
                if covered != full or transversal_ranks(presentation, n)[-1] != r:
                    continue
                check(presentation, n)
        print('General realization through n =', n, ':', count, 'presentations PASS', flush=True)
    for presentation, n in (((7,), 3), ((3, 12), 4), ((3, 4), 3),
                            ((1033, 4612, 2084, 1064, 560, 352, 448, 4114), 13)):
        check(presentation, n, padding=True)
    print('General realization:', count, 'instances,', queries,
          'full-state ranks,', reverse, 'reverse maps PASS', flush=True)


if __name__ == '__main__':
    audit_nontransversal()
    for k in range(4, 13):
        audit_directed_comb(k)
    audit_general_realization()
