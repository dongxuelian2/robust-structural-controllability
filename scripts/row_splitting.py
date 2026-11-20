"""Remove the row-size restriction from exact cotransversal control reductions."""
from bounded_core_solver import matching_number
from hardness_gadgets import cotransversal_graph, protect_coloops


def full_rank_presentation(presentation, n):
    """Compress any transversal presentation to rank-many rows, preserving T.

    Alternating reachability from unmatched rows identifies a row family
    that can match its entire neighborhood. Replace that family by singleton
    rows on this neighborhood. Zero rows are allowed in the input.
    """
    rows = tuple(presentation)
    full = (1 << n) - 1
    if n < 1 or any(row < 0 or row & ~full for row in rows):
        raise ValueError('Use nonnegative n-bit rows.')
    left_mate, right_mate = [-1] * len(rows), [-1] * n
    for start in range(len(rows)):
        queue, seen, predecessor, free = [start], {start}, {}, None
        for left in queue:
            mask = rows[left]
            while mask:
                bit = mask & -mask
                mask ^= bit
                right = bit.bit_length() - 1
                if right in predecessor:
                    continue
                predecessor[right] = left
                other = right_mate[right]
                if other < 0:
                    free = right
                    break
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
            if free is not None:
                break
        if free is None:
            continue
        right = free
        while right >= 0:
            left = predecessor[right]
            old = left_mate[left]
            left_mate[left], right_mate[right] = right, left
            right = old
    reachable_rows = {j for j, mate in enumerate(left_mate) if mate < 0}
    reachable_targets, queue = set(), sorted(reachable_rows)
    for left in queue:
        mask = rows[left]
        while mask:
            bit = mask & -mask
            mask ^= bit
            right = bit.bit_length() - 1
            reachable_targets.add(right)
            other = right_mate[right]
            assert other >= 0  # Otherwise the matching was not maximum.
            if other not in reachable_rows:
                reachable_rows.add(other)
                queue.append(other)
    return (tuple(1 << v for v in sorted(reachable_targets)) +
            tuple(row for j, row in enumerate(rows) if j not in reachable_rows))


def ternary_presentation(presentation, n):
    """Split large rows, returning rows, new ground size and a reversal log."""
    rows = list(full_rank_presentation(presentation, n))
    full = (1 << n) - 1
    if not 0 < len(rows) < n or any(row <= 0 or row & ~full for row in rows):
        raise ValueError('Use nonempty rows on n elements, with 0<r<n.')
    if matching_number(tuple(rows) + (0,) * (n - len(rows)), full) != len(rows):
        raise ValueError('The presentation must have full row rank.')
    covered = 0
    for row in rows:
        covered |= row
    if covered != full:
        raise ValueError('The transversal matroid must have no loops.')
    steps = []
    for start in range(len(rows)):
        i = start
        while rows[i].bit_count() > 3:
            old = rows[i]
            first = old & -old
            second = (old ^ first) & -(old ^ first)
            part = first | second
            steps.append((i, old))
            rows[i] = part | (1 << n)
            rows.append((old ^ part) | (1 << n))
            i = len(rows) - 1
            n += 1
    return tuple(rows), n, tuple(steps)


def dual_rank_oracle(presentation, n):
    """Full-row-rank presentation; return a cached dual-matroid rank oracle."""
    rows = tuple(presentation) + (0,) * (n - len(presentation))
    r, full, cache = len(presentation), (1 << n) - 1, {}

    def rank(selected):
        mask = sum(1 << e for e in selected)
        if mask not in cache:
            cache[mask] = mask.bit_count() + matching_number(rows, full ^ mask) - r
        return cache[mask]

    return rank


def avoid_presentation_element(presentation, n, selected, x):
    """Avoid x at no extra cost when its column has at least two row neighbors.

    The selected set must robustly span the dual transversal matroid. The
    output spans robustly in its deletion of x; that deletion need not itself
    be cotransversal. This operation cannot in general be composed to avoid
    a prescribed set of columns.
    """
    presentation = tuple(presentation)
    if not 0 <= x < n or sum(bool(row >> x & 1) for row in presentation) < 2:
        raise ValueError('The column must have at least two row neighbors.')
    rows = full_rank_presentation(presentation, n)
    rank = dual_rank_oracle(rows, n)
    d, selected = n - len(rows), set(selected)

    def robust(support):
        return rank(support) == d and all(rank(support - {e}) == d for e in support)

    if not selected <= set(range(n)) or not robust(selected):
        raise ValueError('The selected set must be robust.')
    if x not in selected:
        return frozenset(selected)
    selected.remove(x)
    if robust(selected):
        return frozenset(selected)
    replacement = next((w for w in range(n) if w != x and w not in selected
                        and robust(selected | {w})), None)
    if replacement is None:
        raise AssertionError('The tight-Hall-set deletion exchange failed.')
    return frozenset(selected | {replacement})


def normalize_row_splits(presentation, n, selected):
    """Map a robust set after all row splits to a no-larger original robust set."""
    rows, size, steps = ternary_presentation(presentation, n)
    selected = set(selected)
    d = n - len(full_rank_presentation(presentation, n))

    def robust(rank, support):
        return rank(support) == d and all(rank(support - {e}) == d for e in support)

    if not selected <= set(range(size)) or not robust(dual_rank_oracle(rows, size), selected):
        raise ValueError('The selected split elements must be robust.')
    for i, old in reversed(steps):
        rows = list(rows[:-1])
        rows[i] = old
        size -= 1
        if size not in selected:
            continue
        selected.remove(size)
        rank = dual_rank_oracle(rows, size)
        if robust(rank, selected):
            continue
        replacement = next((w for w in range(size) if w not in selected and
                            robust(rank, selected | {w})), None)
        if replacement is None:
            raise AssertionError('The row-splitting exchange failed.')
        selected.add(replacement)
    return frozenset(selected)


def cotransversal_reduction(presentation, n):
    """Exact subcubic realization with arbitrary original row sizes."""
    rows, size, _ = ternary_presentation(presentation, n)
    return cotransversal_graph(rows, size)


def cotransversal_from_placement(presentation, n, instance, selected):
    """Reverse a robust control placement, also after odd-subdivision padding."""
    rows, size, _ = ternary_presentation(presentation, n)
    counts = {}
    for v in set(selected):
        e = instance.projection[v]
        if e >= 0:
            counts[e] = counts.get(e, 0) + 1
    support = protect_coloops(dual_rank_oracle(rows, size), range(size), counts)
    if len(support) > sum(counts.values()):
        raise ValueError('The placement is not robust after parallel projection.')
    return normalize_row_splits(presentation, n, support)
