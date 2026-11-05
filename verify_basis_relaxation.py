"""Exact certificates for the fundamental-assignment relaxation obstruction."""
from collections import Counter
from itertools import combinations
from random import Random
from assignment_transfer import assignment_rank, robust_assignment, normalize_dominated
from bounded_core_solver import matching_number
from verify_gammoid_obstructions import dual_ranks, robust, transversal_ranks


def assignment_ranks(n, options):
    rows = [sum(1 << e for e, option in enumerate(options) if option >> j & 1)
            for j in range(n)]
    return transversal_ranks(rows, len(options))


def fundamental_model(k, t, actual, basis):
    d = t + k
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
    return expected, exchanges


def family(k):
    t, d = 3 * k + 1, 4 * k + 1
    original = tuple(((1 << k) - 1) | (1 << (k + i))
                     for i in range(t) for _ in range(2))
    actual = original + tuple(1 << j for j in range(d))
    basis = set(range(2 * k)) | {2 * i for i in range(k, t)}
    expected, exchanges = fundamental_model(k, t, actual, basis)
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


def sparse_family(k):
    t, d = 3 * k, 4 * k
    original = []
    for group in range(t):
        u, v = group % k, (group + 1) % k
        original.extend(((1 << u) | (1 << (k + group)),
                         (1 << v) | (1 << (k + group))))
    original = tuple(original)
    assert len(set(original)) == 6 * k
    assert all(option.bit_count() == 2 for option in original)
    degrees = [sum(option >> j & 1 for option in original) for j in range(d)]
    assert degrees == [6] * k + [2] * (3 * k)
    actual = original + tuple(1 << j for j in range(d))
    basis = set(range(2 * k)) | {2 * i for i in range(k, t)}
    expected, exchanges = fundamental_model(k, t, actual, basis)
    selected = set(range(2 * k, 2 * t)) | {2 * t}
    assert len(selected) == d + 1
    assert robust_assignment(d, expected, selected)
    assert assignment_rank(d, actual, selected) == 3 * k
    assert robust_assignment(d, original, range(2 * t))
    # A path partition attaining the degree bound pc(H)>=|V(H)|-2k.
    paths = [(2 * k + i, i, 3 * k + i) for i in range(k)]
    paths += [(k + i,) for i in range(k)]
    assert len(paths) == 2 * k
    assert sorted(v for path in paths for v in path) == list(range(d))
    for path in paths:
        assert all((1 << u) | (1 << v) in original for u, v in zip(path, path[1:]))
    if k == 3:
        full = set(range(len(actual)))
        for u, v in combinations(range(len(actual)), 2):
            assert assignment_rank(d, actual, {u, v}) == 2
            assert assignment_rank(d, actual, full - {u, v}) == d
    print('Arity-two family k =', k, ': true OPT', 6 * k, '; relaxed OPT', 4 * k + 1,
          ';', exchanges, 'fundamental exchanges PASS', flush=True)


def prototype():
    presentation = (1033, 4612, 2084, 1064, 560, 352, 448, 4114)
    n, ranks = 13, dual_ranks(presentation, 13)
    d, counts = ranks[-1], Counter()
    robust_masks = [mask for mask in range(1 << n) if robust(ranks, mask)]
    assert min(mask.bit_count() for mask in robust_masks) == 7
    assert min(mask.bit_count() for mask in robust_masks if mask >> 5 & 1) == 8
    for part in combinations(range(n), d):
        basis = sum(1 << v for v in part)
        if ranks[basis] != d:
            continue
        rows = [(1 << b) | sum(1 << e for e in range(n) if not basis >> e & 1
                               and ranks[(basis ^ (1 << b)) | (1 << e)] == d)
                for b in part]
        relaxation = transversal_ranks(rows, n)
        value = min(mask.bit_count() for mask in range(1 << n) if robust(relaxation, mask))
        assert value == (7 if basis >> 5 & 1 else 6)
        counts[value] += 1
        if basis == 79:
            assert value == 6 and robust(relaxation, 215)
            assert ranks[215 ^ 1] == 4
    assert counts == {6: 396, 7: 117}
    print('Nontransversal prototype: all 513 bases;', dict(counts), 'PASS', flush=True)


def degree_basis_family(k):
    t, d = 3 * k + 1, 4 * k + 1
    actual = tuple(((1 << k) - 1) | (1 << (k + i))
                   for i in range(t) for _ in range(2))
    actual += tuple(1 << j for j in range(d))
    base_rows = []
    for i in range(t):
        row = (3 << (2 * i)) | (1 << (2 * t + k + i))
        row |= ((1 << k) - 1) << (2 * t)
        base_rows.extend((row, row))
        members = {e for e in range(len(actual)) if row >> e & 1}
        assert assignment_rank(d, actual, members) == k + 1
        assert all(assignment_rank(d, actual, members | {e}) == k + 2
                   for e in range(len(actual)) if e not in members)
    if k == 1:
        assert dual_ranks(base_rows, len(actual)) == assignment_ranks(d, actual)
    bad = set(range(2 * k)) | {2 * i for i in range(k, t)}
    copies = 2 * t + 1
    origin = list(range(len(actual)))
    for e in sorted(bad):
        origin.extend([e] * (copies - 1))
    classes = [sum(1 << v for v, original in enumerate(origin) if original == e)
               for e in range(len(actual))]
    rows = [sum(classes[e] for e in range(len(actual)) if row >> e & 1)
            for row in base_rows]
    for e in sorted(bad):
        rows.extend([classes[e]] * (copies - 1))
    n, full = len(origin), (1 << len(origin)) - 1

    def rank(selected):
        return assignment_rank(d, actual, {origin[e] for e in selected})

    weights = [sum(row >> e & 1 for row in rows) for e in range(n)]
    assert all(weights[e] == copies + 1 if origin[e] in bad else weights[e] <= 2 * t
               for e in range(n))
    basis = set()
    for e in sorted(range(n), key=lambda e: (-weights[e], e)):
        if rank(basis | {e}) > len(basis):
            basis.add(e)
    assert len(basis) == d and {origin[e] for e in basis} == bad
    ordered, relaxed = sorted(basis), []
    for e in range(n):
        relaxed.append((1 << ordered.index(e)) if e in basis else
                       sum(1 << j for j, b in enumerate(ordered)
                           if rank((basis - {b}) | {e}) == d))
    witness = set(range(2 * k, 2 * t))
    assert robust_assignment(d, relaxed, witness)
    original_model, _ = fundamental_model(k, t, actual, bad)
    assert all(relaxed[e] == original_model[origin[e]] for e in range(n))
    padded_rows = tuple(rows) + (0,) * d
    assert matching_number(padded_rows, full) == len(rows)
    rng = Random(20261105 + k)
    supports = [set(), basis, witness, set(range(n)), set(range(n)) - basis]
    supports += [{e for e in range(n) if rng.randrange(2)} for _ in range(20)]
    for selected in supports:
        complement = full ^ sum(1 << e for e in selected)
        observed = len(selected) + matching_number(padded_rows, complement) - len(rows)
        assert observed == rank(selected)
    print('Maximum-degree basis k =', k, ':', n, 'elements;', len(rows),
          'dual rows; true OPT', 2 * t, '; relaxed OPT', d + 1, 'PASS', flush=True)


if __name__ == '__main__':
    for k in range(1, 13):
        family(k)
    for k in range(3, 13):
        sparse_family(k)
    prototype()
    for k in range(1, 4):
        degree_basis_family(k)
