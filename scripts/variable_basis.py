"""Private-basis / bicircular family for investigating hardness without forcing."""
from collections import Counter
from hardness_gadgets import Instance, protect_coloops


def private_realize(n, options):
    """General option neighborhoods as nonempty n-bit masks, with privates."""
    options = tuple(options)
    if n < 1 or any(option <= 0 or option >> n for option in options):
        raise ValueError('Use nonempty option neighborhoods on range(n).')
    covered = 0
    for option in options:
        covered |= option
    if covered != (1 << n) - 1:
        raise ValueError('Every job must have an option.')
    m = len(options)
    rows = [0] * (n + 2 * m)
    for j, option in enumerate(options):
        c = n + m + j
        for t in [n + j] + [v for v in range(n) if option >> v & 1]:
            rows[c] |= 1 << t
            rows[t] |= 1 << c
    return tuple(rows)


def realize(n, edges):
    """State b_v=v; edge resource t_e=n+e; center c_e=n+m+e.

    Use a simple graph on range(n), without isolated vertices. Connected
    source graphs give strongly connected state graphs. No state self-loops.
    """
    edges = tuple(sorted({tuple(sorted(e)) for e in edges}))
    if n < 2 or any(u == v or u < 0 or v >= n for u, v in edges):
        raise ValueError('Use a simple graph on range(n), n>=2.')
    if {v for e in edges for v in e} != set(range(n)):
        raise ValueError('The source graph has an isolated vertex.')
    m = len(edges)
    rows = [0] * (n + 2 * m)
    for j, (u, v) in enumerate(edges):
        c = n + m + j
        for t in (u, v, n + j):
            rows[c] |= 1 << t
            rows[t] |= 1 << c
    return tuple(rows), edges


def sparse_realize(n, edges):
    """Split high-degree job states into parallel copies; maximum degree three."""
    original, edges = realize(n, edges)
    m = len(edges)
    rows = [0] * len(original)
    labels = ([('b', v) for v in range(n)] + [('t', j) for j in range(m)] +
              [('c', j) for j in range(m)])
    centers = set(range(n + m, n + 2 * m))
    projection = list(range(n + m)) + [-1] * m

    def link(u, v):
        rows[u] |= 1 << v
        rows[v] |= 1 << u

    for j in range(m):
        link(n + j, n + m + j)
    for v in range(n):
        incident = [j for j, e in enumerate(edges) if v in e]
        previous = v
        for i, j in enumerate(incident):
            if i:
                c, t = len(rows), len(rows) + 1
                rows.extend((0, 0))
                labels.extend((('split_center', v, i), ('b_copy', v, i)))
                projection.extend((-1, v))
                centers.add(c)
                link(previous, c)
                link(c, t)
                previous = t
            link(previous, n + m + j)
    return Instance(tuple(rows), tuple(labels), edges, frozenset(centers),
                    (), (), {}, {}, tuple(projection))


def improve_cover(n, edges, cover=None):
    """Merge paths/cycles and close paths until no listed move applies.

    Each component is (is_cycle, ordered_vertices). Lexicographically
    (#components, #path_components) decreases; the procedure is polynomial.
    """
    components = [(False, [v]) for v in range(n)] if cover is None else [
        (cycle, list(vertices)) for cycle, vertices in cover]
    edges = tuple(edges)
    while True:
        location = {v: j for j, (_, part) in enumerate(components) for v in part}
        changed = False
        for u, v in edges:
            i, j = location[u], location[v]
            ci, a = components[i]
            cj, b = components[j]
            eu, ev = u in (a[0], a[-1]), v in (b[0], b[-1])
            if i == j:
                if not ci and len(a) >= 3 and {u, v} == {a[0], a[-1]}:
                    components[i] = True, a
                    changed = True
                    break
                continue
            merged = None
            if not ci and not cj and eu and ev:
                left = a if a[-1] == u else a[::-1]
                right = b if b[0] == v else b[::-1]
                merged = left + right
            elif not ci and cj and eu:
                left = a if a[-1] == u else a[::-1]
                cut = b.index(v)
                merged = left + b[cut:] + b[:cut]
            elif ci and not cj and ev:
                cut = a.index(u)
                right = b if b[0] == v else b[::-1]
                merged = a[cut + 1:] + a[:cut + 1] + right
            elif ci and cj:
                cut_a, cut_b = a.index(u), b.index(v)
                merged = (a[cut_a + 1:] + a[:cut_a + 1] +
                          b[cut_b:] + b[:cut_b])
            if merged is not None:
                components = [part for k, part in enumerate(components)
                              if k not in (i, j)] + [(False, merged)]
                changed = True
                break
        if not changed:
            return tuple((cycle, tuple(part)) for cycle, part in components)


def cover_from_placement(n, edges, selected):
    """For a robust placement return <= |S|-n paths/cycles, using chosen edges."""
    edges = tuple(sorted({tuple(sorted(e)) for e in edges}))
    selected = set(selected)
    chosen = [edge for j, edge in enumerate(edges) if n + j in selected]
    result = improve_cover(n, chosen)
    assert len(result) <= len(selected & set(range(n + len(edges)))) - n
    return result


def cover_from_parallel(n, edges, instance, selected):
    """Extract a cover from parallel copies, including protected coloops."""
    counts = Counter(instance.projection[v] for v in set(selected)
                     if instance.projection[v] >= 0)
    rows, edges = realize(n, edges)
    from bounded_core_solver import matching_number
    full = (1 << len(rows)) - 1
    nu = matching_number(rows, full)
    cache = {}

    def rank(resources):
        mask = sum(1 << v for v in resources)
        if mask not in cache:
            cache[mask] = len(resources) + matching_number(rows, full ^ mask) - nu
        return cache[mask]

    normalized = protect_coloops(rank, range(n + len(edges)), counts)
    assert len(normalized) <= sum(counts.values())
    result = cover_from_placement(n, edges, normalized)
    return result


def placement_from_cover(n, edges, cover):
    """Lift a path/cycle cover to at most n+#components actuators."""
    edges = tuple(sorted({tuple(sorted(e)) for e in edges}))
    cover = improve_cover(n, edges, cover)
    edge_index = {edge: n + j for j, edge in enumerate(edges)}
    selected, isolated = set(), set()
    for cycle, part in cover:
        if len(part) == 1:
            isolated.add(part[0])
            continue
        for u, v in zip(part, part[1:]):
            selected.add(edge_index[tuple(sorted((u, v)))])
        if cycle:
            selected.add(edge_index[tuple(sorted((part[-1], part[0])))])
            selected.add(part[0])
        else:
            selected.update((part[0], part[-1]))
    for v in isolated:
        edge = next(edge for edge in edges if v in edge and
                    any(w not in isolated for w in edge))
        selected.update((v, edge_index[edge]))
    assert len(selected) <= n + len(cover)
    return frozenset(selected)
