"""Exact certificates for general source-SCC and cotransversal universality."""
from itertools import combinations, combinations_with_replacement, product
from random import Random
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from quota_cotransversal import quota_cotransversal_graph, quota_from_placement, sparse_control_graph
from row_splitting import dual_rank_oracle
from verify_gammoid_obstructions import dual_ranks, robust
from verify_scc_cover import girth


def input_first_matchable(rows, selected):
    """Direct augmented matching, reserving dedicated input rows first.

    This avoids repeatedly rerouting a nearly complete input matching on
    the large padded fixtures. It uses no dual-rank or projection formula.
    """
    size = len(rows)
    selected = sorted(selected)
    augmented = [1 << v for v in selected] + list(rows)
    mate = [-1] * size
    for i, v in enumerate(selected):
        mate[v] = i

    def augment(u, seen):
        remaining = augmented[u]
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            v = bit.bit_length() - 1
            if v in seen:
                continue
            seen.add(v)
            if mate[v] < 0 or augment(mate[v], seen):
                mate[v] = u
                return True
        return False

    matched = len(selected)
    if matched == size:
        return True
    for u in range(len(selected), len(augmented)):
        matched += augment(u, set())
        if matched == size:
            return True
    return False


def check(presentation, n, groups, rng, length=1, exhaustive_source=True):
    instance = quota_cotransversal_graph(presentation, n, groups, length)
    rows, size = instance.rows, len(instance.rows)
    full, d = (1 << size) - 1, n - len(instance.presentation)
    assert matching_number(rows, full) == size - d
    sources = source_components(rows)
    assert set(sources) == {sum(1 << v for v in cycle) for cycle in instance.source_cycles}
    assert all(block.bit_count() >= 4 * length for block in sources)
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1) for v in range(size)]
    single = [mask for mask in incoming if mask.bit_count() == 1]
    assert all(incoming) and len(single) == len(set(single))
    assert max(mask.bit_count() for mask in incoming) <= 3
    assert max(row.bit_count() for row in rows) <= 3
    neighbors = [{v for v in range(size) if rows[u] >> v & 1 or rows[v] >> u & 1}
                 for u in range(size)]
    color, queue = {0: 0}, [0]
    for u in queue:
        assert 2 <= len(neighbors[u]) <= 3 and u not in neighbors[u]
        for v in neighbors[u]:
            if v not in color:
                color[v] = 1 - color[u]
                queue.append(v)
            assert color[v] != color[u]
    assert len(color) == size
    assert size == instance.core_vertices + (length - 1) * instance.core_edges + 4 * length * instance.leaves
    q = sum(row.bit_count() for row in presentation)
    assert size <= 13 * q * length
    ranks = dual_ranks(instance.presentation, n)
    reps = instance.representatives
    rank_checks = forward = reverse = failures = 0

    def actual_rank(selected):
        return len(selected) + matching_number(rows, full ^ sum(1 << v for v in selected)) - (size - d)

    def actual_robust(selected):
        nonlocal failures
        for v in selected:
            failures += 1
            survivors = selected - {v}
            actual = input_first_matchable(rows, survivors)
            if size <= 40:
                assert actual == augmented_matchable(rows, survivors)
            if not actual:
                return False
        return bool(selected)

    original_masks = (range(1 << n) if exhaustive_source else
                      sorted({0, (1 << n) - 1, 223} | {rng.randrange(1 << n) for _ in range(32)}))
    optimum = n + 1
    for mask in original_masks:
        chosen = {e for e in range(n) if mask >> e & 1}
        lifted = {reps[e] for e in chosen}
        assert actual_rank(lifted) == ranks[mask]
        rank_checks += 1
        if robust(ranks, mask) and all(len(chosen & set(group)) >= 2 for group in groups):
            assert actual_robust(lifted)
            forward += 1
            optimum = min(optimum, len(chosen))

    exhaustive = size <= 12
    if exhaustive:
        supports = [{v for v in range(size) if mask >> v & 1} for mask in range(1 << size)]
    else:
        supports = [{v for v in range(size) if rng.random() < probability}
                    for probability in (.3, .65, .9) for _ in range(2 if length > 1 else 4)]
    if length == 1 and size <= 40:
        supports += [set(range(size)) - {v} for v in range(size)]
    else:
        supports += [set(range(size)) - {v} for v in rng.sample(range(size), min(4, size))]
    target_optimum = size + 1
    for selected in supports:
        if any(len(selected & {v for v in range(size) if block >> v & 1}) < 2 for block in sources):
            continue
        valid = actual_robust(selected)
        if len(selected) == size - 1:
            assert valid
        if not valid:
            continue
        result = quota_from_placement(instance, selected)
        assert {e for e in range(n) if reps[e] in selected} <= result
        assert len(result) <= len(selected)
        assert robust(ranks, sum(1 << e for e in result))
        assert all(len(result & set(group)) >= 2 for group in groups)
        target_optimum = min(target_optimum, len(selected))
        reverse += 1
    if exhaustive:
        assert target_optimum == optimum

    split_rank = dual_rank_oracle(instance.split_rows, instance.split_size)
    for _ in range(10):
        selected = set(rng.sample(range(size), rng.randrange(min(size, n + 3) + 1)))
        projected = {instance.projection[v] for v in selected if instance.projection[v] >= 0}
        assert actual_rank(selected) == split_rank(projected)
        rank_checks += 1
    if length > 1:
        base = quota_cotransversal_graph(presentation, n, groups)
        assert girth(rows) == length * girth(base.rows)
    return len(supports), forward, reverse, rank_checks, failures


def control_fixtures():
    """Positive-deficiency sources with several SCCs and ungrouped states."""
    rng = Random(202611052)
    fixtures = [(18, 33, 24, 36, 0, 0), (66, 132, 193, 80, 160, 200, 0, 0)]
    for n in (4, 5, 6):
        accepted = 0
        while accepted < 4:
            rows = tuple(sum(1 << v for v in range(n) if rng.randrange(3) == 0)
                         for _ in range(n))
            if any(block.bit_count() < 2 for block in source_components(rows)):
                continue
            fixtures.append(rows)
            accepted += 1
    totals, count = [0] * 5, 0
    for rows in fixtures:
        n = len(rows)
        mapped = sparse_control_graph(rows)
        groups = mapped.groups
        source_ranks = dual_rank_oracle(mapped.presentation, n)
        d = n - len(mapped.presentation)
        # Independently compare the source's actual input-failure model.
        for mask in range(1 << n):
            selected = {v for v in range(n) if mask >> v & 1}
            if not all(len(selected & set(group)) >= 2 for group in groups):
                continue
            source_valid = all(augmented_matchable(rows, selected - {v}) for v in selected)
            expected = source_ranks(selected) == d and all(
                source_ranks(selected - {v}) == d for v in selected)
            assert source_valid == expected
        for length in ((1, 3) if n >= 7 else (1,)):
            counts = check(rows, n, groups, rng, length)
            totals = [a + b for a, b in zip(totals, counts)]
            count += 1
    print('Additional directed control fixtures:', count, 'instances;', totals,
          'placements / forward / reverse / ranks / failures PASS', flush=True)
    return count, totals


if __name__ == '__main__':
    rng = Random(202611051)
    totals = [0] * 5
    instances = 0

    def run(presentation, n, groups, length=1, exhaustive_source=True):
        global instances
        counts = check(presentation, n, groups, rng, length, exhaustive_source)
        for i, count in enumerate(counts):
            totals[i] += count
        instances += 1

    for n in (2, 3):
        groups = [(group,) for size in range(2, n + 1) for group in combinations(range(n), size)]
        for k in range(1, 4):
            for presentation in combinations_with_replacement(range(1 << n), k):
                if any(not any(row >> e & 1 for row in presentation) for e in range(n)):
                    continue
                for partition in groups:
                    run(presentation, n, partition)
        print('General quota presentations through', n, 'elements:', instances, 'instances;',
              totals, 'placements / forward / reverse / ranks / failures PASS', flush=True)

    control_sources = 0
    for n in (2, 3):
        for rows in product(range(1 << n), repeat=n):
            source = source_components(rows)
            if any(block.bit_count() < 2 for block in source):
                continue
            groups = tuple(tuple(v for v in range(n) if block >> v & 1) for block in source)
            mapped = sparse_control_graph(rows)
            assert mapped.groups == groups
            run(rows, n, groups)
            control_sources += 1
    print('All feasible directed sources through three states:', control_sources, 'PASS', flush=True)

    fixtures = [
        ((2, 1, 8, 4), 4, ((0, 1), (2, 3))),
        ((2, 1, 7, 12), 4, ((0, 1),)),
        ((15,), 4, ((0, 1),)),
        ((1, 2, 12), 4, ((0, 1),)),
        ((63, 7, 56), 6, ((0, 1, 2),)),
        ((31,) * 5, 5, ((0, 1),)),
    ]
    for presentation, n, groups in fixtures:
        for length in (1, 3, 9):
            run(presentation, n, groups, length)
            print('Quota fixture:', n, 'elements; length', length, 'PASS', flush=True)
    prototype = (1033, 4612, 2084, 1064, 560, 352, 448, 4114)
    run(prototype, 13, ((0, 1, 2),), exhaustive_source=False)
    run(prototype, 13, ((0, 1, 2),), 3, exhaustive_source=False)
    print('Certified nontransversal prototype with a proper quota subset PASS', flush=True)
    print('TOTAL:', instances, 'instances;', totals,
          'placements / forward / reverse / ranks / failures PASS', flush=True)
    control_fixtures()
