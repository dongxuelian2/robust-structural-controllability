"""Direct state/input matching certificates for the source-SCC cover reduction."""
from itertools import combinations, combinations_with_replacement
from collections import Counter, deque
from random import Random
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from scc_cover import cover_graph, extract_covers, lift_covers, complete_basis, pad_cover


def minimum_cover(n, sets):
    for size in range(1, len(sets) + 1):
        for chosen in combinations(range(len(sets)), size):
            union = 0
            for i in chosen:
                union |= sets[i]
            if union == (1 << n) - 1:
                return frozenset(chosen)
    raise ValueError('The source is not coverable.')


def check(instance, exhaustive, rng):
    rows, n = instance.rows, instance.jobs
    size, full = len(rows), (1 << len(rows)) - 1
    nu = matching_number(rows, full)
    assert size - nu == 2 * n
    sources = source_components(rows)
    expected_sources = [{j, n + j} for j in range(n)]
    expected_sources += [set(instance.options[0]) | set(instance.options[1])]
    assert {frozenset(v for v in range(size) if block >> v & 1) for block in sources} == {
        frozenset(block) for block in expected_sources}
    assert augmented_matchable(rows, instance.basis)
    # A column with one incoming neighbor has a different such neighbor
    # from every other one-neighbor column; other columns have degree >=2.
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1)
                for v in range(size)]
    single = [mask for mask in incoming if mask.bit_count() == 1]
    assert all(incoming) and len(single) == len(set(single))
    extras = sorted(set(range(size)) - instance.basis)
    if exhaustive:
        supports = [instance.basis | {v for i, v in enumerate(extras) if mask >> i & 1}
                    for mask in range(1 << len(extras))]
    else:
        supports = [instance.basis | {v for v in extras if rng.randrange(2)}
                    for _ in range(30)]
    optimum, valid_count, queries = size + 1, 0, 0
    for selected in supports:
        reached = 0
        for v in selected - instance.basis:
            reached |= instance.neighborhoods[v]
        quota = sum(v in selected for group in instance.options for v in group) >= 2
        predicted = quota and reached == (1 << (2 * n)) - 1
        actual = quota
        for failed in instance.basis:
            queries += 1
            if not augmented_matchable(rows, selected - {failed}):
                actual = False
                break
        assert actual == predicted, (instance.sets, selected)
        if not actual:
            continue
        first, second = extract_covers(instance, selected)
        recovered = lift_covers(instance, first, second)
        assert len(recovered) <= len(selected)
        assert all(augmented_matchable(rows, recovered - {v}) for v in recovered)
        valid_count += 1
        optimum = min(optimum, len(selected))
    cover = minimum_cover(n, instance.sets)
    witness = lift_covers(instance, cover)
    assert all(augmented_matchable(rows, witness - {v}) for v in witness)
    if exhaustive:
        assert optimum == 2 * (n + len(cover))
    return len(supports), valid_count, queries


def girth(rows):
    size = len(rows)
    neighbors = [{v for v in range(size) if row >> v & 1 or rows[v] >> u & 1}
                 for u, row in enumerate(rows)]
    best = size + 1
    for start in range(size):
        distance, parent, queue = {start: 0}, {start: -1}, deque([start])
        while queue:
            u = queue.popleft()
            for v in neighbors[u]:
                if v not in distance:
                    distance[v], parent[v] = distance[u] + 1, u
                    queue.append(v)
                elif parent[u] != v and parent[v] != u:
                    best = min(best, distance[u] + distance[v] + 1)
    return best


def check_optional(n, sets, rng):
    raw, instance = cover_graph(n, sets, optional=False), cover_graph(n, sets)
    rows, size = instance.rows, len(instance.rows)
    # This exact row-multiset identity certifies all ranks, not just samples:
    # the new columns are free transversal summands, hence dual loops.
    assert Counter(r for r in rows if r) == Counter(r for r in raw.rows if r) + Counter(
        1 << v for v in range(len(raw.rows), size))
    sources = source_components(rows)
    assert set(sources) == {sum(1 << v for v in group) for group in instance.source_cycles}
    assert all(group.bit_count() >= 4 for group in sources)
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1)
                for v in range(size)]
    single = [mask for mask in incoming if mask.bit_count() == 1]
    assert all(incoming) and len(single) == len(set(single))
    # Simple T and source sizes >=3 prove that the full ground minus any
    # one element is robust. Check those placements directly in small cases.
    supports = [set(range(size)) - {v} for v in range(size)]
    if size <= 12:
        supports += [{v for v in range(size) if mask >> v & 1} for mask in range(1 << size)]
    else:
        supports += [{v for v in range(size) if rng.random() < probability}
                     for probability in (.35, .6, .9) for _ in range(20)]
    valid_count, optimum = 0, size
    for selected in supports:
        mask = sum(1 << v for v in selected)
        valid = all((mask & group).bit_count() >= 2 for group in sources) and all(
            augmented_matchable(rows, selected - {v}) for v in selected)
        if len(selected) == size - 1:
            assert valid
        if not valid:
            continue
        normalized = complete_basis(instance, selected)
        assert len(normalized) <= len(selected)
        assert all(augmented_matchable(rows, normalized - {v}) for v in normalized)
        first, second = extract_covers(instance, selected)
        recovered = lift_covers(instance, first, second)
        assert len(recovered) <= len(selected)
        assert all(augmented_matchable(rows, recovered - {v}) for v in recovered)
        valid_count += 1
        optimum = min(optimum, len(selected))
    if size <= 12:
        assert optimum == 2 * (n + len(minimum_cover(n, sets)))
    print('Optional source cycles:', n, 'jobs;', size, 'states;', len(supports),
          'placements;', valid_count, 'reverse maps PASS', flush=True)
    return raw, instance


def check_padding(raw, base, length, rng):
    instance = pad_cover(base, length)
    rows, size = instance.rows, len(instance.rows)
    full, deficiency = (1 << size) - 1, 2 * instance.jobs
    nu = size - deficiency
    edges = sum(bool(base.rows[u] >> v & 1 or base.rows[v] >> u & 1)
                for u in range(len(base.rows)) for v in range(u + 1, len(base.rows)))
    assert size == len(base.rows) + (length - 1) * edges
    assert girth(rows) == length * girth(base.rows)
    assert matching_number(rows, full) == nu
    assert set(source_components(rows)) == {
        sum(1 << v for v in group) for group in instance.source_cycles}
    color, queue = {0: 0}, [0]
    for u in queue:
        neighbors = {v for v in range(size) if rows[u] >> v & 1 or rows[v] >> u & 1}
        assert 2 <= len(neighbors) <= 3
        for v in neighbors:
            if v not in color:
                color[v] = 1 - color[u]
                queue.append(v)
            assert color[v] != color[u]
    assert len(color) == size
    assert max(row.bit_count() for row in rows) <= 3
    assert max(sum(row >> v & 1 for row in rows) for v in range(size)) <= 3
    assert all(not row >> v & 1 for v, row in enumerate(rows))
    witness = lift_covers(instance, minimum_cover(instance.jobs, instance.sets))
    basis = set(instance.basis)
    supports = [set(), basis, witness, set(range(size))]
    supports += [set(rng.sample(range(size), rng.randrange(1, deficiency + 3)))
                 for _ in range(15)]
    supports += [(basis - {rng.choice(sorted(basis))}) | {rng.randrange(size)}
                 for _ in range(15)]
    for selected in supports:
        projected = {instance.projection[v] for v in selected if instance.projection[v] >= 0}
        observed = len(selected) + matching_number(
            rows, full ^ sum(1 << v for v in selected)) - nu
        old_full = (1 << len(raw.rows)) - 1
        expected = len(projected) + matching_number(
            raw.rows, old_full ^ sum(1 << v for v in projected)) - (len(raw.rows) - deficiency)
        assert observed == expected
    assert all(augmented_matchable(rows, witness - {v}) for v in witness)
    # The old full set minus any basis element stays robust under padding.
    # Its source loop selections exercise completion without a full basis.
    for omitted in sorted(basis)[:4]:
        selected = set(range(len(base.rows))) - {omitted}
        selected |= {v for v in range(len(base.rows), size) if rng.random() < .1}
        first, second = extract_covers(instance, selected)
        recovered = lift_covers(instance, first, second)
        assert len(recovered) <= len(selected)
        assert all(augmented_matchable(rows, recovered - {v}) for v in recovered)
    print('Odd padding length', length, ':', size, 'states; girth', girth(rows),
          ';', len(supports), 'full-state dual ranks PASS', flush=True)


if __name__ == '__main__':
    rng = Random(20261105)
    instances = placements = reverse = queries = 0
    for n in (1, 2, 3):
        for q in range(1, 4):
            for sets in combinations_with_replacement(range(1, 1 << n), q):
                if any(not any(s >> j & 1 for s in sets) for j in range(n)):
                    continue
                for sparse in (False, True):
                    instance = cover_graph(n, sets, sparse, optional=False)
                    exhaustive = len(instance.rows) - 2 * n <= 12
                    a, b, c = check(instance, exhaustive, rng)
                    instances += 1
                    placements += a
                    reverse += b
                    queries += c
        print('SCC cover through', n, 'jobs:', instances, 'instances;', placements,
              'placements;', reverse, 'reverse maps;', queries, 'state/input tests PASS', flush=True)
    fixtures = [tuple(combinations(range(5), 2)),
                ((0, 1), (1, 2), (2, 3), (3, 0)),
                ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))]
    for edges in fixtures:
        vertices = sorted({v for edge in edges for v in edge})
        sets = tuple(sum(1 << j for j, edge in enumerate(edges) if v in edge)
                     for v in vertices)
        instance = cover_graph(len(edges), sets, optional=False)
        rows, size = instance.rows, len(instance.rows)
        degree = [sum(bool(rows[u] >> v & 1 or rows[v] >> u & 1)
                      for v in range(size)) for u in range(size)]
        assert set(degree) <= {2, 3}
        assert max(row.bit_count() for row in rows) <= 3
        assert max(sum(row >> v & 1 for row in rows) for v in range(size)) <= 3
        assert all(not row >> v & 1 for v, row in enumerate(rows))
        a, b, c = check(instance, False, rng)
        cover = minimum_cover(len(edges), sets)
        print('Subcubic fixture:', len(edges), 'source edges;', size, 'states; OPT',
              2 * (len(edges) + len(cover)), ';', c, 'matching tests PASS', flush=True)
        if len(edges) == 10:
            completed = cover_graph(len(edges), sets)
            completed_edges = sum(bool(completed.rows[u] >> v & 1 or completed.rows[v] >> u & 1)
                                  for u in range(len(completed.rows))
                                  for v in range(u + 1, len(completed.rows)))
            assert len(completed.rows) == 14 * len(edges) - 2 * len(vertices)
            assert completed_edges == 18 * len(edges) - 2 * len(vertices)
            check_padding(instance, completed, 3, rng)
    for n, sets in ((1, (1,)), (1, (1, 1)), (2, (1, 2)), (3, (3, 6, 5)),
                    (4, (3, 5, 8, 6))):
        raw, base = check_optional(n, sets, rng)
        for length in (3, 9):
            check_padding(raw, base, length, rng)
