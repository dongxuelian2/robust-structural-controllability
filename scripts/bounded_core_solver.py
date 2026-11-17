"""Exact single-actuator-failure solver, exponential only for unbounded d.

Input: rows[i] is the bitmask of out-neighbors of state i, 0 <= i < n.
Output: an optimum tuple of distinct actuator states, or None if infeasible.
This reference implementation uses no external packages and does not mutate input.
"""
from itertools import combinations


def matching_number(rows, allowed):
    n = len(rows)
    mate = [-1] * n
    count = 0
    for start in range(n):
        queue = [start]
        seen_left = {start}
        predecessor = {}
        free = None
        for left in queue:
            mask = rows[left] & allowed
            while mask:
                bit = mask & -mask
                mask ^= bit
                right = bit.bit_length() - 1
                if right in predecessor:
                    continue
                predecessor[right] = left
                if mate[right] == -1:
                    free = right
                    break
                other = mate[right]
                if other not in seen_left:
                    seen_left.add(other)
                    queue.append(other)
            if free is not None:
                break
        if free is None:
            continue
        # Follow the alternating path backward via old matched right endpoints.
        left_mate = {left: right for right, left in enumerate(mate) if left >= 0}
        right = free
        while True:
            left = predecessor[right]
            old = left_mate.get(left)
            mate[right] = left
            if old is None:
                break
            right = old
        count += 1
    return count


def source_components(rows):
    n = len(rows)
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    remaining = (1 << n) - 1
    result = []
    while remaining:
        i = (remaining & -remaining).bit_length() - 1
        component = sum(1 << j for j in range(n)
                        if reach[i] & (1 << j) and reach[j] & (1 << i))
        remaining &= ~component
        if not any(row & component for j, row in enumerate(rows)
                   if not component & (1 << j)):
            result.append(component)
    return result


def solve(rows):
    rows = tuple(rows)
    n = len(rows)
    if n == 0:
        raise ValueError('The state graph must be nonempty.')
    full = (1 << n) - 1
    if any(not isinstance(row, int) or row < 0 or row & ~full for row in rows):
        raise ValueError('Each row must be a nonnegative n-bit integer.')
    sources = source_components(rows)
    if any(c.bit_count() < 2 for c in sources):
        return None
    nu = matching_number(rows, full)
    d = n - nu
    rank_cache = {}
    def rank(s):
        if s not in rank_cache:
            rank_cache[s] = s.bit_count() + matching_number(rows, full ^ s) - nu
        return rank_cache[s]
    best = full
    for size in range(min(n, 2*d) + 1):
        if size > best.bit_count():
            break
        for indices in combinations(range(n), size):
            core = sum(1 << i for i in indices)
            if rank(core) != d or any(rank(core ^ (1 << i)) != d for i in indices):
                continue
            candidate = core
            for component in sources:
                need = max(0, 2 - (core & component).bit_count())
                available = component & ~core
                for _ in range(need):
                    bit = available & -available
                    candidate |= bit
                    available ^= bit
            if candidate.bit_count() < best.bit_count():
                best = candidate
    return tuple(i for i in range(n) if best & (1 << i))
