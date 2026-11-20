"""Exact weighted spanning by flat saturation and bounded search."""
from functools import cache
from itertools import combinations
from math import comb

from weighted_fault_core import weighted_input


def obstruction_bound(rank, failures, size_limit, packed_bases=1):
    """Upper bound on the hyperplanes to avoid in a same-status exchange."""
    general = comb(rank + failures - 1, rank - 1)
    if packed_bases == 1:
        return general
    missing = min(rank, (failures + 1) // packed_bases)
    packed = packed_bases * sum(
        comb(rank, omitted) * comb(size_limit - rank, omitted - 1)
        for omitted in range(1, missing + 1)
    )
    return 1 + min(general, packed)


def flat_core(elements, weights, rank, ambient_rank, size_limit, obstructions):
    """Retain cheap nonloops of one status using rank-indexed flat caps."""
    thresholds = [0]
    for _ in range(ambient_rank):
        thresholds.append(obstructions * thresholds[-1] + size_limit)

    retained = []
    discarded = 0
    ordered = sorted(elements, key=lambda element: (weights[element], element))

    def in_saturated_proper_flat(element):
        for dimension in range(1, ambient_rank):
            threshold = thresholds[dimension]
            if len(retained) < threshold:
                break
            for subset in combinations(retained, dimension):
                basis = frozenset(subset)
                if rank(basis) != dimension or rank(basis | {element}) != dimension:
                    continue
                if sum(rank(basis | {e}) == dimension for e in retained) >= threshold:
                    return True
        return False

    for index, element in enumerate(ordered):
        if len(retained) == thresholds[ambient_rank]:
            discarded += len(ordered) - index
            break
        if in_saturated_proper_flat(element):
            discarded += 1
        else:
            retained.append(element)
    return frozenset(retained), discarded


def _packs_bases(candidate, rank, ambient_rank, count):
    """Test disjoint base packing by disjoint-subset recursion."""
    @cache
    def packs(available, remaining):
        if len(available) < remaining * ambient_rank:
            return False
        if remaining == 1:
            return rank(available) == ambient_rank
        return any(
            packs(available - basis, remaining - 1)
            for subset in combinations(sorted(available), ambient_rank)
            if rank(basis := frozenset(subset)) == ambient_rank
        )

    return packs(candidate, count)


def solve(weights, rank_oracle, k, vulnerable=None, *, size_limit=None, packed_bases=1):
    """Return cost, selected set, core and counts for the requested restrictions.

    The optional cap and packing requirement apply to the selected set.
    An unrestricted call minimizes cost over all resilient sets.
    """
    weights, ground, vulnerable = weighted_input(weights, k, vulnerable)
    if size_limit is not None and (type(size_limit) is not int or size_limit < 0):
        raise ValueError("size_limit must be a nonnegative integer")
    if type(packed_bases) is not int or not 1 <= packed_bases <= k + 1:
        raise ValueError("packed_bases must be between 1 and k+1")

    @cache
    def small_rank(elements):
        return rank_oracle(elements)

    ambient_rank = rank_oracle(ground)
    limit = min(len(ground), (k + 1) * ambient_rank)
    if size_limit is not None:
        limit = min(limit, size_limit)
    empty = frozenset()
    stats = dict(discarded=0, candidates=0, cached_queries=0)
    if ambient_rank == 0:
        return 0, empty, empty, stats
    if limit < ambient_rank * packed_bases:
        return None, None, empty, stats

    def ordered(elements):
        return sorted(elements, key=lambda e: (weights[e], e))

    nonloops = frozenset(e for e in ground if small_rank(frozenset({e})))
    safe = nonloops - vulnerable
    exposed = nonloops & vulnerable

    def rank(elements):
        if len(elements) <= ambient_rank:
            return small_rank(elements)
        basis = empty
        for element in sorted(elements):
            extension = basis | {element}
            if small_rank(extension) == len(extension):
                basis = extension
                if len(basis) == ambient_rank:
                    break
        return len(basis)

    def result(witness, core):
        stats["cached_queries"] = small_rank.cache_info().misses
        cost = None if witness is None else sum(weights[e] for e in witness)
        return cost, witness, core, stats

    if k == 0:
        basis = empty
        for element in ordered(nonloops):
            extension = basis | {element}
            if small_rank(extension) == len(extension):
                basis = extension
                if len(basis) == ambient_rank:
                    return result(basis, basis)

    if ambient_rank == 1:
        core = frozenset(ordered(safe)[:limit] + ordered(exposed)[:limit])
        candidates = []
        if safe and len(core) >= packed_bases:
            cheapest_safe = ordered(safe)[0]
            companions = ordered(core - {cheapest_safe})[:packed_bases - 1]
            candidates.append(frozenset([cheapest_safe, *companions]))
        if limit >= k + 1 and len(exposed) >= k + 1:
            candidates.append(frozenset(ordered(exposed)[:k + 1]))
        witness = min(candidates, key=lambda selected: (
            sum(weights[e] for e in selected), len(selected), tuple(sorted(selected))
        ), default=None)
        return result(witness, core)

    bound = obstruction_bound(ambient_rank, k, limit, packed_bases)
    safe_core, safe_discards = flat_core(safe, weights, small_rank, ambient_rank, limit, bound)
    exposed_core, exposed_discards = flat_core(
        exposed, weights, small_rank, ambient_rank, limit, bound
    )
    core = safe_core | exposed_core
    stats["discarded"] = safe_discards + exposed_discards
    best = None
    witness = None
    core_order = sorted(core)
    for size in range(ambient_rank * packed_bases, min(limit, len(core)) + 1):
        for elements in combinations(core_order, size):
            stats["candidates"] += 1
            cost = sum(weights[e] for e in elements)
            if best is not None and cost >= best:
                continue
            candidate = frozenset(elements)
            failing = sorted(candidate & vulnerable)
            if any(rank(candidate - frozenset(failed)) < ambient_rank
                   for failed in combinations(failing, min(k, len(failing)))):
                continue
            if packed_bases > 1 and not _packs_bases(candidate, rank, ambient_rank, packed_bases):
                continue
            best, witness = cost, candidate
    return result(witness, core)


def solve_independence(weights, independent, k, vulnerable=None, *,
                       size_limit=None, packed_bases=1):
    """Solve with a Boolean independence oracle on frozensets of element IDs."""
    def rank(elements):
        basis = frozenset()
        for element in sorted(elements):
            extension = basis | {element}
            if independent(extension):
                basis = extension
        return len(basis)

    return solve(weights, rank, k, vulnerable,
                 size_limit=size_limit, packed_bases=packed_bases)
