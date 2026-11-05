"""Direct state/input matching certificates for the source-SCC cover reduction."""
from itertools import combinations, combinations_with_replacement
from random import Random
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from scc_cover import cover_graph, extract_covers, lift_covers


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
    assert augmented_matchable(rows, instance.forced)
    # A column with one incoming neighbor has a different such neighbor
    # from every other one-neighbor column; other columns have degree >=2.
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1)
                for v in range(size)]
    single = [mask for mask in incoming if mask.bit_count() == 1]
    assert all(incoming) and len(single) == len(set(single))
    extras = sorted(set(range(size)) - instance.forced)
    if exhaustive:
        supports = [instance.forced | {v for i, v in enumerate(extras) if mask >> i & 1}
                    for mask in range(1 << len(extras))]
    else:
        supports = [instance.forced | {v for v in extras if rng.randrange(2)}
                    for _ in range(30)]
    optimum, valid_count, queries = size + 1, 0, 0
    for selected in supports:
        reached = 0
        for v in selected - instance.forced:
            reached |= instance.neighborhoods[v]
        quota = sum(v in selected for group in instance.options for v in group) >= 2
        predicted = quota and reached == (1 << (2 * n)) - 1
        actual = quota
        for failed in instance.forced:
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


if __name__ == '__main__':
    rng = Random(20261105)
    instances = placements = reverse = queries = 0
    for n in (1, 2, 3):
        for q in range(1, 4):
            for sets in combinations_with_replacement(range(1, 1 << n), q):
                if any(not any(s >> j & 1 for s in sets) for j in range(n)):
                    continue
                for sparse in (False, True):
                    instance = cover_graph(n, sets, sparse)
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
        instance = cover_graph(len(edges), sets)
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
