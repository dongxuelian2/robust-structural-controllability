"""Exact checks of dominated-resource normalization and the branching transfer."""
from itertools import combinations_with_replacement
from random import Random
from assignment_transfer import (assignment_from_placement, assignment_graph,
                                 branching_presentation, normalize_dominated,
                                 robust_assignment)
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from hardness_gadgets import subdivide_ports
from verify_gammoid_obstructions import dual_ranks, robust, transversal_ranks


def resource_ranks(n, options):
    rows = [sum(1 << e for e, option in enumerate(options) if option >> j & 1)
            for j in range(n)]
    return transversal_ranks(rows, len(options))


def audit_exchange():
    count = exchanges = 0
    for n in range(1, 4):
        for m in range(n + 1, 7):
            for options in combinations_with_replacement(range(1, 1 << n), m):
                ranks = resource_ranks(n, options)
                if ranks[-1] != n or not robust(ranks, (1 << m) - 1):
                    continue
                original_robust = [robust(ranks, mask) for mask in range(1 << m)]
                for added in range(1, 1 << n):
                    if not any(added & base == added for base in options):
                        continue
                    extended = resource_ranks(n, options + (added,))
                    for mask in range(1 << m, 1 << (m + 1)):
                        if not robust(extended, mask):
                            continue
                        removed = mask ^ (1 << m)
                        assert original_robust[removed] or any(
                            not removed >> w & 1 and original_robust[removed | (1 << w)]
                            for w in range(m))
                        exchanges += 1
                count += 1
        print('Dominated-resource exchanges through', n, 'jobs:', count,
              'sources,', exchanges, 'exchanges PASS', flush=True)
    return count, exchanges


def audit_branching():
    rng = Random(51629)
    count = rank_queries = reverse = 0

    def check(n, options, exhaustive=False, pad=False):
        nonlocal count, rank_queries, reverse
        options = tuple(options)
        m = len(options)
        original = resource_ranks(n, options)
        if original[-1] != n or not robust(original, (1 << m) - 1):
            return
        rows, relaxed = branching_presentation(n, options)
        size = len(relaxed)
        assert size - len(rows) == n and max(row.bit_count() for row in rows) <= 3
        assert size - n <= sum(option.bit_count() for option in options)
        ranks = dual_ranks(rows, size) if exhaustive else None
        upper = resource_ranks(n, relaxed) if exhaustive else None
        full = (1 << size) - 1

        def rank(mask):
            return (ranks[mask] if ranks is not None else mask.bit_count() +
                    matching_number(rows + (0,) * n, full ^ mask) - len(rows))

        source_masks = (range(1 << m) if exhaustive else
                        [(1 << m) - 1] + [rng.randrange(1 << m) for _ in range(30)])
        for mask in source_masks:
            assert rank(mask) == original[mask]
        supports = list(range(1 << size)) if exhaustive else []
        supports += [sum(1 << v for v in range(size) if rng.randrange(2))
                     for _ in range(30)]
        for mask in supports:
            rank_queries += 1
            if exhaustive:
                assert rank(mask) <= upper[mask]
            if rank(mask) != n or any(rank(mask ^ (1 << v)) != n
                                      for v in range(size) if mask >> v & 1):
                continue
            selected = {v for v in range(size) if mask >> v & 1}
            normalized = normalize_dominated(n, options, relaxed, selected)
            assert len(normalized) <= len(selected)
            assert robust(original, sum(1 << v for v in normalized))
            reverse += 1

        graph = assignment_graph(n, options)
        if pad:
            graph = subdivide_ports(graph)
        state_count = len(graph.rows)
        assert max(row.bit_count() for row in graph.rows) <= 3
        assert all(not row >> v & 1 for v, row in enumerate(graph.rows))
        assert len(source_components(graph.rows)) == 1
        assert len(graph.rows) - matching_number(graph.rows, (1 << state_count) - 1) == n
        for mask in source_masks:
            if robust(original, mask):
                selected = {e for e in range(m) if mask >> e & 1}
                assert all(augmented_matchable(graph.rows, selected - {e})
                           for e in selected)
        for _ in range(8):
            selected = {v for v in range(state_count) if rng.randrange(2)}
            if not selected or not all(augmented_matchable(graph.rows, selected - {v})
                                       for v in selected):
                continue
            normalized = assignment_from_placement(n, options, graph, selected)
            assert len(normalized) <= len(selected)
            assert robust(original, sum(1 << v for v in normalized))
            reverse += 1
        count += 1

    for n in range(1, 4):
        for m in range(n + 1, 5):
            for options in combinations_with_replacement(range(1, 1 << n), m):
                check(n, options, exhaustive=True)
        print('Branching transfer through', n, 'jobs:', count, 'sources PASS', flush=True)
    # A source whose optimum is six, strictly above rank plus one.
    options = (3, 3, 3, 12, 12, 12)
    ranks = resource_ranks(4, options)
    assert min(mask.bit_count() for mask in range(1 << 6) if robust(ranks, mask)) == 6
    check(4, options, exhaustive=True, pad=True)
    # Large original neighborhoods actually require splitting row centers.
    check(5, (31, 31, 1, 2, 4, 8, 16))
    check(6, (63, 63, 1, 2, 4, 8, 16, 32), pad=True)
    for n in range(4, 8):
        for _ in range(12):
            # Two full resources plus singleton resources guarantee feasibility.
            options = [(1 << n) - 1, (1 << n) - 1] + [1 << j for j in range(n)]
            options += [rng.randrange(1, 1 << n) for _ in range(2)]
            check(n, options)
    print('Branching transfer:', count, 'sources,', rank_queries, 'network supports,',
          reverse, 'reverse maps PASS', flush=True)


if __name__ == '__main__':
    audit_exchange()
    audit_branching()
