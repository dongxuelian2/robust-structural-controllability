"""Independent exact check of the private-basis / hyperforest formula."""
from itertools import combinations_with_replacement
from functools import reduce
from random import Random
from hardness_gadgets import protect_coloops
from variable_basis import private_realize
from audit_hardness import augmented_matchable


def assignment_ranks(n, options):
    """Enumerate partial injections from jobs to available resources."""
    job_rows = [1 << v for v in range(n)]
    for j, option in enumerate(options):
        for v in range(n):
            if option >> v & 1:
                job_rows[v] |= 1 << (n + j)
    reachable = {0}
    for row in job_rows:
        for mask in tuple(reachable):
            remaining = row & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                reachable.add(mask | bit)
    ranks = [mask.bit_count() if mask in reachable else 0
             for mask in range(1 << (n + len(options)))]
    for j in range(n + len(options)):
        for mask in range(1 << (n + len(options))):
            if mask >> j & 1:
                ranks[mask] = max(ranks[mask], ranks[mask ^ (1 << j)])
    return ranks


def hyperforest_rank(n, options, chosen):
    """Rank from the strong Hall inequalities on every job subset."""
    neighbors = [sum(1 << j for j, option in enumerate(options)
                     if chosen >> j & 1 and option >> v & 1)
                 for v in range(n)]
    unions = [0] * (1 << n)
    independent = [True] * (1 << n)
    best = 0
    for mask in range(1, 1 << n):
        bit = mask & -mask
        unions[mask] = unions[mask ^ bit] | neighbors[bit.bit_length() - 1]
        independent[mask] = unions[mask].bit_count() > mask.bit_count()
        remaining = mask
        while remaining and independent[mask]:
            bit = remaining & -remaining
            remaining ^= bit
            independent[mask] = independent[mask ^ bit]
        if independent[mask]:
            best = max(best, mask.bit_count())
    return best


def partitions(items):
    """All set partitions, represented by disjoint nonempty bit masks."""
    if not items:
        yield ()
        return
    bit = items & -items
    for rest in partitions(items ^ bit):
        yield (bit,) + rest
        for j in range(len(rest)):
            yield rest[:j] + (rest[j] | bit,) + rest[j + 1:]


def partition_deficit(n, options, chosen):
    neighbors = [sum(1 << j for j, option in enumerate(options)
                     if chosen >> j & 1 and option >> v & 1)
                 for v in range(n)]
    return max(len(parts) - sum(sum(bool(edge & part) for part in parts) >= 2
                                for edge in neighbors)
               for parts in partitions(chosen))


def check(n, options):
    full_jobs = (1 << n) - 1
    if not options or reduce(int.__or__, options) != full_jobs:
        return False
    ranks = assignment_ranks(n, options)
    rows = private_realize(n, options) if n <= 3 else None
    actual = n + len(options)
    for selected, rank in enumerate(ranks):
        if rank != n:
            continue
        support = {v for v in range(n + len(options)) if selected >> v & 1}
        critical = sum(ranks[selected ^ (1 << v)] < n for v in support)
        normalized = protect_coloops(
            lambda states: ranks[sum(1 << v for v in states)],
            range(n + len(options)), support)
        normalized_mask = sum(1 << v for v in normalized)
        assert len(normalized) <= len(support) + critical
        assert all(ranks[normalized_mask ^ (1 << v)] == n for v in normalized)
        remaining = selected
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            if ranks[selected ^ bit] != n:
                break
        else:
            actual = min(actual, selected.bit_count())
        if rows is not None:
            robust = all(ranks[selected ^ (1 << v)] == n for v in support)
            assert augmented_matchable(rows, support)
            assert robust == all(augmented_matchable(rows, support - {v})
                                 for v in support)
    expected = n + len(options)
    for chosen in range(1, 1 << len(options)):
        covered = 0
        for j, option in enumerate(options):
            if chosen >> j & 1:
                covered |= option
        if covered != full_jobs:
            continue
        rank = hyperforest_rank(n, options, chosen)
        deficit = chosen.bit_count() - rank
        assert deficit == partition_deficit(n, options, chosen), (n, options, chosen)
        expected = min(expected, n + deficit)
    assert actual == expected, (n, options, actual, expected)
    return True


if __name__ == '__main__':
    count = 0
    for n in range(1, 5):
        choices = [mask for mask in range(1, 1 << n) if mask.bit_count() <= 3]
        for m in range(1, 5):
            for options in combinations_with_replacement(choices, m):
                count += check(n, options)
        print('Exhaustive private-basis instances through n =', n, ':', count)
    rng = Random(20261104)
    random_count = 0
    for n in (5, 6):
        for _ in range(50):
            options = tuple(sum(1 << v for v in rng.sample(range(n), rng.randrange(1, 4)))
                            for _ in range(rng.randrange(3, 7)))
            random_count += check(n, options)
    print('Random instances:', random_count)
    print('PASS: direct robust matching = n + min_C (|C| - hyperforest rank).')
