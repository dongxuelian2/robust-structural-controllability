"""Unit-cost uniform robust assignments as subcubic control instances."""
from bounded_core_solver import matching_number
from hardness_gadgets import cotransversal_graph, protect_coloops


def assignment_rank(n, options, selected):
    """Maximum assignment cardinality using the selected resource indices."""
    size = max(n, len(options))
    rows = [sum(1 << e for e, option in enumerate(options) if option >> j & 1)
            for j in range(n)] + [0] * (size - n)
    return matching_number(rows, sum(1 << e for e in selected))


def robust_assignment(n, options, selected):
    selected = set(selected)
    return bool(selected) and assignment_rank(n, options, selected) == n and all(
        assignment_rank(n, options, selected - {e}) == n for e in selected)


def normalize_dominated(n, original, extended, selected):
    """Remove dominated added resources without increasing selected cardinality.

    All resources have unit cost and are vulnerable. Original resources must
    already form a feasible robust instance. Each added neighborhood must be
    contained in an original neighborhood; identical neighborhoods are allowed.
    """
    original, extended = tuple(original), tuple(extended)
    m = len(original)
    if (n < 1 or extended[:m] != original or
            any(option <= 0 or option >> n for option in extended) or
            any(not any(option & base == option for base in original)
                for option in extended[m:])):
        raise ValueError('Added resources must be dominated by original resources.')
    selected = set(selected)
    if (not selected <= set(range(len(extended))) or
            not robust_assignment(n, original, range(m)) or
            not robust_assignment(n, extended, selected)):
        raise ValueError('Both the original instance and selected set must be robust.')
    for e in sorted(selected - set(range(m))):
        removed = selected - {e}
        if robust_assignment(n, extended, removed):
            selected = removed
            continue
        replacement = next((w for w in range(m) if w not in selected and
                            robust_assignment(n, extended, removed | {w})), None)
        if replacement is None:
            raise AssertionError('The tight-Hall-set exchange failed.')
        selected = removed | {replacement}
    return frozenset(selected)


def branching_presentation(n, options):
    """Dual presentation of disjoint resource trees sharing their job sinks.

    Elements 0,...,m-1 are original resource roots, m,...,m+n-1 are jobs.
    Later elements are branch vertices. Return rows and their relaxed
    assignment neighborhoods (reachable jobs), in element-index order.
    """
    options = tuple(options)
    if n < 1 or any(option <= 0 or option >> n for option in options):
        raise ValueError('Use nonempty resource neighborhoods on range(n).')
    if not robust_assignment(n, options, range(len(options))):
        raise ValueError('The original uniform robust assignment must be feasible.')
    m = len(options)
    relaxed = list(options) + [1 << j for j in range(n)]
    rows = []
    for e, option in enumerate(options):
        jobs = [j for j in range(n) if option >> j & 1]
        current = e
        while len(jobs) > 2:
            branch = len(relaxed)
            relaxed.append(sum(1 << j for j in jobs[1:]))
            rows.append((1 << current) | (1 << (m + jobs[0])) | (1 << branch))
            current, jobs = branch, jobs[1:]
        rows.append((1 << current) | sum(1 << (m + j) for j in jobs))
    return tuple(rows), tuple(relaxed)


def assignment_graph(n, options):
    """Return an exact, strongly connected subcubic control realization."""
    rows, relaxed = branching_presentation(n, options)
    return cotransversal_graph(rows, len(relaxed))


def assignment_from_placement(n, options, instance, selected):
    """Recover a no-larger robust original assignment from a robust placement.

    `instance` may also have been padded with `subdivide_ports`; its projection
    refers to the elements of `branching_presentation`.
    """
    rows, relaxed = branching_presentation(n, options)
    size, r = len(relaxed), len(rows)
    full = (1 << size) - 1
    counts = {}
    for v in set(selected):
        e = instance.projection[v]
        if e >= 0:
            counts[e] = counts.get(e, 0) + 1

    def rank(elements):
        mask = sum(1 << e for e in elements)
        return len(elements) + matching_number(rows + (0,) * n, full ^ mask) - r

    support = protect_coloops(rank, range(size), counts)
    if len(support) > sum(counts.values()):
        raise ValueError('The placement is not robust after parallel projection.')
    return normalize_dominated(n, options, relaxed, support)
