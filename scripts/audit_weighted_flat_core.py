"""Finite checks for weighted flat saturation and disjoint base packing."""
from itertools import combinations
from math import comb
from random import Random

from audit_weighted_core import feasible_masks, make_rank, members
from audit_weighted_independent import binary_ranks
from weighted_flat_core import flat_core, obstruction_bound, solve, solve_independence


def check_instance(weights, rank_mask, k, vulnerable_mask):
    n = len(weights)
    ground = (1 << n) - 1
    r = rank_mask(ground)
    vulnerable = members(n, vulnerable_mask)
    feasible = feasible_masks(n, rank_mask, k, vulnerable_mask)
    bases = [sum(1 << e for e in subset) for subset in combinations(range(n), r)
             if rank_mask(sum(1 << e for e in subset)) == r]
    checks = 0

    def rank(elements):
        return rank_mask(sum(1 << e for e in elements))

    for count in range(1, k + 2):
        packed_unions = set()
        for group in combinations(bases, count):
            union = 0
            for basis in group:
                if union & basis:
                    break
                union |= basis
            else:
                packed_unions.add(union)
        if r == 0:
            packed_unions = {0}
        for limit in {n, max(0, (n + r) // 2), r * count}:
            expected = {
                candidate: sum(weights[e] for e in members(n, candidate))
                for candidate in feasible
                if candidate.bit_count() <= limit
                and any(union & candidate == union for union in packed_unions)
            }
            actual, witness, core, _ = solve(
                iter(weights), rank, k, iter(vulnerable),
                size_limit=limit, packed_bases=count,
            )
            assert actual == min(expected.values(), default=None)
            if witness is not None:
                mask = sum(1 << e for e in witness)
                assert expected[mask] == actual and witness <= core
            else:
                assert actual is None
            checks += 1
    return checks


def check_discard_witnesses():
    labels = [0] * 7 + [1] * 2 + list(range(2, 18))
    weights = tuple(7 if e == 8 else e for e in range(len(labels)))
    _, rank = make_rank('classes', len(labels), labels=labels, trunc=2)
    ground = frozenset(range(len(labels)))
    core, discarded = flat_core(ground, weights, rank, 2, 4, 2)
    assert len(core) == 12 and discarded == len(labels) - 12
    assert core & frozenset(range(7)) == frozenset(range(4))
    assert solve(weights, rank, 1)[:2] == (15, frozenset({0, 1, 7, 8}))


def check_packed_bound():
    r, k, count = 5, 6, 7
    assert obstruction_bound(r, k, count * r, count) == 1 + count * r
    assert 1 + count * r < comb(r + k - 1, k)

    r, k, count = 4, 4, 5
    labels = list(range(r)) * count
    ground = frozenset(range(len(labels)))
    _, partition_rank = make_rank('classes', len(labels), labels=labels, trunc=r)
    assert obstruction_bound(r, k, len(labels), count) == 21 < comb(r + k - 1, k)
    assert solve((1,) * len(labels), partition_rank, k,
                 packed_bases=count)[:2] == (len(labels), ground)

    # The full ground has two disjoint bases, but its cheapest resilient
    # subset uses three classes once each and cannot pack the second base.
    labels = [0, 0, 1, 1, 2]
    weights = (1, 8, 1, 8, 1)
    _, rank = make_rank('classes', len(labels), labels=labels, trunc=2)
    assert solve(weights, rank, 1)[0] == 3
    assert solve(weights, rank, 1, packed_bases=2)[0] == 11

    def independent(elements):
        assert isinstance(elements, frozenset)
        return rank(elements) == len(elements)

    assert solve_independence(weights, independent, 1, packed_bases=2)[0] == 11


def check_query_cache():
    labels = list(range(16))
    _, rank = make_rank('classes', len(labels), labels=labels, trunc=2)
    calls = set()
    ground = frozenset(range(len(labels)))

    def once(elements):
        if elements != ground:
            assert len(elements) <= 2 and elements not in calls
            calls.add(elements)
        return rank(elements)

    cost, _, _, stats = solve((0,) * len(labels), once, 2)
    assert cost == 0 and stats['cached_queries'] == len(calls)


def main():
    rng = Random(20261120)
    checks = 0
    for n in range(9):
        for _ in range(18):
            ranks = binary_ranks([rng.randrange(8) for _ in range(n)])
            weights = tuple(rng.choice([0, 1, 2, 6, 100, 10**40]) for _ in range(n))
            checks += check_instance(weights, ranks.__getitem__, rng.randrange(3),
                                     rng.randrange(1 << n))
    vamos, _ = make_rank('sparse_paving', 8, r=4, forbidden=[
        (0, 1, 2, 3), (0, 1, 4, 5), (0, 1, 6, 7), (2, 3, 4, 5), (2, 3, 6, 7),
    ])
    for k in range(3):
        for vulnerable in (0, 0b10110101, 255):
            checks += check_instance(tuple(range(8)), vamos, k, vulnerable)
    check_discard_witnesses()
    check_packed_bound()
    check_query_cache()
    print(f'PASS: {checks} exact comparisons, including the Vamos matroid.')
    print('PASS: proper and ambient saturation, disjoint base packing, oracle caching.')


if __name__ == '__main__':
    main()
