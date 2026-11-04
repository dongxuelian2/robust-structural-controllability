"""Exact certificates for the fundamental-assignment relaxation obstruction."""
from collections import Counter
from itertools import combinations
from assignment_transfer import assignment_rank, robust_assignment, normalize_dominated
from verify_gammoid_obstructions import dual_ranks, robust, transversal_ranks


def assignment_ranks(n, options):
    rows = [sum(1 << e for e, option in enumerate(options) if option >> j & 1)
            for j in range(n)]
    return transversal_ranks(rows, len(options))


def family(k):
    t, d = 3 * k + 1, 4 * k + 1
    original = tuple(((1 << k) - 1) | (1 << (k + i))
                     for i in range(t) for _ in range(2))
    actual = original + tuple(1 << j for j in range(d))
    basis = set(range(2 * k)) | {2 * i for i in range(k, t)}
    ordered, anchors = sorted(basis), (1 << (2 * k)) - 1
    expected = []
    for e in range(len(actual)):
        if e in basis:
            mask = 1 << ordered.index(e)
        elif e < 2 * t:
            mask = anchors | (1 << (k + e // 2))
        elif e - 2 * t < 2 * k:
            mask = anchors
        else:
            i = e - 2 * t - k
            mask = anchors | (1 << (k + i))
        expected.append(mask)
    assert assignment_rank(d, actual, basis) == d
    exchanges = 0
    for e in range(len(actual)):
        if e in basis:
            continue
        observed = sum(1 << j for j, b in enumerate(ordered)
                       if assignment_rank(d, actual, (basis - {b}) | {e}) == d)
        assert observed == expected[e]
        exchanges += len(basis)
    selected = set(range(2 * k, 2 * t))
    assert len(selected) == d + 1
    assert robust_assignment(d, expected, selected)
    assert assignment_rank(d, actual, selected) == 3 * k + 1
    assert robust_assignment(d, original, range(2 * t))
    # An actual robust set with privates normalizes to all original resources.
    enlarged = set(range(2 * t)) | set(range(2 * t, len(actual), 2))
    assert normalize_dominated(d, original, actual, enlarged) == frozenset(range(2 * t))
    if k == 1:
        rank, relaxation = assignment_ranks(d, actual), assignment_ranks(d, expected)
        assert min(s.bit_count() for s in range(len(rank)) if robust(rank, s)) == 8
        assert min(s.bit_count() for s in range(len(rank)) if robust(relaxation, s)) == 6
        full = len(rank) - 1
        assert all(rank[(1 << u) | (1 << v)] == 2 and
                   rank[full ^ (1 << u) ^ (1 << v)] == d
                   for u, v in combinations(range(len(actual)), 2))
    print('Family k =', k, ': true OPT', 6 * k + 2, '; relaxed OPT', 4 * k + 2,
          ';', exchanges, 'fundamental exchanges PASS', flush=True)


def prototype():
    presentation = (1033, 4612, 2084, 1064, 560, 352, 448, 4114)
    n, ranks = 13, dual_ranks(presentation, 13)
    d, counts = ranks[-1], Counter()
    assert min(mask.bit_count() for mask in range(1 << n) if robust(ranks, mask)) == 7
    for part in combinations(range(n), d):
        basis = sum(1 << v for v in part)
        if ranks[basis] != d:
            continue
        rows = [(1 << b) | sum(1 << e for e in range(n) if not basis >> e & 1
                               and ranks[(basis ^ (1 << b)) | (1 << e)] == d)
                for b in part]
        relaxation = transversal_ranks(rows, n)
        value = min(mask.bit_count() for mask in range(1 << n) if robust(relaxation, mask))
        counts[value] += 1
        if basis == 79:
            assert value == 6 and robust(relaxation, 215)
            assert ranks[215 ^ 1] == 4
    assert counts == {6: 396, 7: 117}
    print('Nontransversal prototype: all 513 bases;', dict(counts), 'PASS', flush=True)


if __name__ == '__main__':
    for k in range(1, 13):
        family(k)
    prototype()
