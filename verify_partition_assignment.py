"""Exact source quotas, full state/input failures and reverse assignment maps."""
from itertools import combinations_with_replacement
from random import Random
from assignment_transfer import assignment_rank, robust_assignment
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from partition_assignment_transfer import (
    partition_assignment_graph, assignment_from_partition_placement)
from verify_private_basis_spares import partitions
from verify_scc_cover import girth


def check(n, options, groups, rng, length=1):
    instance = partition_assignment_graph(n, options, groups, length)
    rows, m = instance.rows, len(options)
    size, full = len(rows), (1 << len(rows)) - 1
    assert matching_number(rows, full) == size - n
    sources = source_components(rows)
    assert set(sources) == {sum(1 << v for v in group) for group in instance.source_cycles}
    assert all(group.bit_count() >= 4 for group in sources)
    color, unseen = {}, set(range(size))
    while unseen:
        start = min(unseen)
        color[start], queue = 0, [start]
        for u in queue:
            neighbors = {v for v in range(size) if rows[u] >> v & 1 or rows[v] >> u & 1}
            assert 2 <= len(neighbors) <= 3 and not rows[u] >> u & 1
            for v in neighbors:
                if v not in color:
                    color[v] = 1 - color[u]
                    queue.append(v)
                assert color[v] != color[u]
        unseen.difference_update(queue)
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1)
                for v in range(size)]
    singleton = [mask for mask in incoming if mask.bit_count() == 1]
    assert all(incoming) and len(set(singleton)) == len(singleton)
    assert max(row.bit_count() for row in rows) <= 3
    assert max(mask.bit_count() for mask in incoming) <= 3
    source_optimum, forward = m + 1, 0
    for mask in range(1 << m):
        chosen = {e for e in range(m) if mask >> e & 1}
        observed = len(chosen) + matching_number(rows, full ^ mask) - (size - n)
        assert observed == assignment_rank(n, options, chosen)
        if (any(len(chosen & set(group)) < 2 for group in groups) or
                not robust_assignment(n, options, chosen)):
            continue
        assert all(augmented_matchable(rows, chosen - {e}) for e in chosen)
        source_optimum = min(source_optimum, len(chosen))
        forward += 1
    assert source_optimum <= m
    exhaustive = size <= 11
    if exhaustive:
        supports = [{v for v in range(size) if mask >> v & 1} for mask in range(1 << size)]
    else:
        supports = [{v for v in range(size) if rng.random() < probability}
                    for probability in (.35, .65, .9) for _ in range(5)]
    if length == 1 and size <= 30:
        supports += [set(range(size)) - {v} for v in range(size)]
    actual_optimum, reverse = size, 0
    for selected in supports:
        mask = sum(1 << v for v in selected)
        valid = all((mask & group).bit_count() >= 2 for group in sources) and all(
            augmented_matchable(rows, selected - {v}) for v in selected)
        if len(selected) == size - 1:
            assert valid
        if not valid:
            continue
        original_selected = selected & set(range(m))
        result = assignment_from_partition_placement(instance, selected)
        assert original_selected <= result
        assert len(result) <= len(selected)
        assert robust_assignment(n, options, result)
        assert all(len(result & set(group)) >= 2 for group in groups)
        assert all(augmented_matchable(rows, result - {e}) for e in result)
        actual_optimum = min(actual_optimum, len(selected))
        reverse += 1
    if exhaustive:
        assert actual_optimum == source_optimum
    if length > 1:
        base = partition_assignment_graph(n, options, groups)
        assert girth(rows) == length * girth(base.rows)
        edges = sum(bool(base.rows[u] >> v & 1 or base.rows[v] >> u & 1)
                    for u in range(len(base.rows)) for v in range(u + 1, len(base.rows)))
        assert size == len(base.rows) + (length - 1) * edges
    else:
        q = sum(option.bit_count() for option in options)
        assert size <= 4 * q - 3 * n
    # Independent full-state rank projection through job splitting and padding.
    base = instance.base_rows
    base_full = (1 << len(base)) - 1
    for _ in range(8):
        selected = set(rng.sample(range(size), rng.randrange(0, min(size, n + 3) + 1)))
        projected = {instance.projection[v] for v in selected if instance.projection[v] >= 0}
        actual = len(selected) + matching_number(rows, full ^ sum(1 << v for v in selected)) - (size - n)
        expected = len(projected) + matching_number(
            base, base_full ^ sum(1 << v for v in projected)) - (len(base) - n)
        assert actual == expected
    return len(supports), forward, reverse


if __name__ == '__main__':
    rng = Random(20261105)
    instances = placements = forward = reverse = 0
    for n in range(1, 4):
        for m in range(2, 5):
            groups = [tuple(tuple(v for v in range(m) if group >> v & 1) for group in partition)
                      for partition in partitions((1 << m) - 1)
                      if all(group.bit_count() >= 2 for group in partition)]
            for options in combinations_with_replacement(range(1, 1 << n), m):
                if not robust_assignment(n, options, range(m)):
                    continue
                for partition in groups:
                    a, b, c = check(n, options, partition, rng)
                    instances += 1
                    placements += a
                    forward += b
                    reverse += c
        print('Partition assignment through', n, 'jobs:', instances, 'instances;',
              placements, 'placements;', forward, 'forward maps;', reverse,
              'reverse maps PASS', flush=True)
    for n, options, groups in (
            (2, (3, 3, 3, 3), ((0, 1), (2, 3))),
            (3, (7, 7, 3, 5, 6), ((0, 1), (2, 3, 4))),
            (4, (15, 15, 7, 11, 13, 14), ((0, 1, 2), (3, 4, 5)))):
        for length in (3, 9):
            a, b, c = check(n, options, groups, rng, length)
            print('Padded partition assignment:', n, 'jobs; length', length,
                  ';', a, 'placements;', b, 'forward;', c, 'reverse PASS', flush=True)
    n = 4
    options = tuple(1 << v for v in range(n)) + tuple(
        (1 << v) | (1 << ((v + 1) % n)) for v in range(n))
    feasible = [mask for mask in range(1 << len(options)) if robust_assignment(
        n, options, {v for v in range(len(options)) if mask >> v & 1})]
    assert min(mask.bit_count() for mask in feasible) == 5
    assert min(mask.bit_count() for mask in feasible if mask & 15 == 15) == 6
    a, b, c = check(n, options, (tuple(range(len(options))),), rng)
    print('Genuine basis choices: cycle OPT five, fixed private basis six;',
          b, 'forward maps;', c, 'reverse maps PASS', flush=True)
