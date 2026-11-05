"""Set Cover via source-SCC quotas, with OPT = 2*(|U| + tau)."""
from dataclasses import dataclass


@dataclass
class CoverInstance:
    rows: tuple
    sets: tuple
    jobs: int
    basis: frozenset
    options: tuple
    neighborhoods: tuple
    owner: dict
    source_cycles: tuple
    projection: tuple


def cover_graph(n, sets, sparse=True, optional=True):
    """Return a loopless bipartite control instance and its cover maps.

    Two copies of the covering instance give OPT=2*(n+tau), without parity
    assumptions or a special case for tau=1. In the sparse version every
    underlying undirected degree is two or three, except jobs with source
    frequency greater than two. Its maximum directed outdegree is three.
    The default adds matroid-loop states to source cycles, preserving OPT
    while making every individual actuator position avoidable.
    """
    sets = tuple(sets)
    full = (1 << n) - 1
    if (n < 1 or not sets or any(s <= 0 or s & ~full for s in sets) or
            any(not any(s >> j & 1 for s in sets) for j in range(n))):
        raise ValueError('Use nonempty subsets covering a nonempty universe.')
    q = len(sets)
    count = 4 * n + 2 * q
    rows, neighborhoods, color = [0] * count, [0] * count, [0] * count
    basis, owner = frozenset(range(2 * n)), {}
    options = tuple(tuple(2 * n + 2 * i + c for i in range(q)) for c in (0, 1))

    def job(c, j):
        return 2 * n + 2 * q + c * n + j

    for c in (0, 1):
        for j in range(n):
            b, mate, target = c * n + j, (1 - c) * n + j, job(c, j)
            rows[mate] |= (1 << b) | (1 << target)
            neighborhoods[b] = neighborhoods[target] = 1 << (c * n + j)
            color[b] = color[target] = c
            owner[target] = c, next(i for i, s in enumerate(sets) if s >> j & 1)

    def attach(parent, c, members, original):
        if len(members) == 1:
            rows[parent] |= 1 << job(c, members[0])
            return
        x, y = len(rows), len(rows) + 1
        rows.extend((1 << y, 1 << x))
        neighborhoods.extend((sum(1 << (c * n + j) for j in members), 0))
        color.extend((c, 1 - c))
        owner[x] = c, original
        rows[parent] |= 1 << x
        middle = len(members) // 2
        attach(y, c, members[:middle], original)
        attach(y, c, members[middle:], original)

    cycle = tuple(range(2 * n, 2 * n + 2 * q))
    for index, v in enumerate(cycle):
        c, original = index % 2, index // 2
        parent = cycle[index - 1]
        rows[parent] |= 1 << v
        color[v] = c
        neighborhoods[v] = sets[original] << (c * n)
        owner[v] = c, original
        members = [j for j in range(n) if sets[original] >> j & 1]
        if sparse:
            attach(parent, c, members, original)
        else:
            rows[parent] |= sum(1 << job(c, j) for j in members)
    source_cycles = [(j, n + j) for j in range(n)] + [cycle]
    projection = list(range(len(rows)))
    if optional:
        for group in range(len(source_cycles)):
            if len(source_cycles[group]) >= 3:
                continue
            old = source_cycles[group]
            u, v = old[-1], old[0]
            x, y = len(rows), len(rows) + 1
            # Relocate the old u-row, including its tree port, to y.
            # The two new columns have singleton rows, so they are T-free.
            rows.extend((1 << y, rows[u]))
            rows[u] = 1 << x
            neighborhoods.extend((0, 0))
            color.extend((1 - color[u], color[u]))
            projection.extend((-1, -1))
            source_cycles[group] = old + (x, y)
    assert all(color[u] != color[v] for u, row in enumerate(rows)
               for v in range(len(rows)) if row >> v & 1)
    return CoverInstance(tuple(rows), sets, n, basis, options,
                         tuple(neighborhoods), owner, tuple(source_cycles),
                         tuple(projection))


def pad_cover(instance, length):
    """Oddly subdivide the underlying graph, preserving the covering OPT.

    Source-cycle ports move with their rows; new source states are M-loops.
    A bidirected internal gate edge gets a bidirected path. Other arcs get
    a directed path whose internal vertices occur in directed two-cycles.
    These latter operations add parallel copies and M-loops. All original
    source cycles must have at least three vertices. Apply once to a
    completed construction; choose the desired total subdivision length.
    """
    rows, cycles, neighborhoods, projection, owner = pad_resource_network(
        instance.rows, instance.source_cycles, instance.neighborhoods,
        instance.projection, instance.owner, length)
    return CoverInstance(rows, instance.sets, instance.jobs, instance.basis,
                         instance.options, neighborhoods, owner, cycles, projection)


def pad_resource_network(original, source_cycles, original_neighborhoods,
                         original_projection, original_owner, length):
    """Pad source cycles and private gate trees; return their exact metadata.

    Every bidirected edge outside source cycles must be an internal gate
    with one endpoint's only out-neighbor its partner. This shared local
    construction also applies to partition-constrained assignment networks.
    """
    if length < 1 or length % 2 == 0 or any(len(c) < 3 for c in source_cycles):
        raise ValueError('Use a positive odd length and source cycles of size >=3.')
    if length == 1:
        return (tuple(original), tuple(source_cycles), tuple(original_neighborhoods),
                tuple(original_projection), dict(original_owner))
    rows = [0] * len(original)
    neighborhoods, projection = list(original_neighborhoods), list(original_projection)
    owner, row_image, cycle_arcs, cycles = dict(original_owner), {}, set(), []

    def state(origin=None):
        v = len(rows)
        rows.append(0)
        neighborhoods.append(0 if origin is None else original_neighborhoods[origin])
        projection.append(-1 if origin is None else original_projection[origin])
        if origin in original_owner:
            owner[v] = original_owner[origin]
        return v

    for cycle in source_cycles:
        expanded = []
        for i, u in enumerate(cycle):
            v = cycle[(i + 1) % len(cycle)]
            internal = [state() for _ in range(length - 1)]
            path = [u] + internal + [v]
            for a, b in zip(path, path[1:]):
                rows[a] |= 1 << b
            row_image[u] = internal[-1]
            cycle_arcs.add((u, v))
            expanded.extend(path[:-1])
        cycles.append(tuple(expanded))
    handled = set(cycle_arcs)
    for u, row in enumerate(original):
        for v in range(len(original)):
            if not row >> v & 1 or (u, v) in handled:
                continue
            if original[v] >> u & 1:
                # One endpoint of an internal gate has only its partner
                # as out-neighbor. Even path positions copy that endpoint.
                x, y = (u, v) if original[u] == 1 << v else (v, u)
                assert original[x] == 1 << y
                internal = [state(x if i % 2 == 0 else None)
                            for i in range(1, length)]
                path = [x] + internal + [y]
                for a, b in zip(path, path[1:]):
                    rows[a] |= 1 << b
                    rows[b] |= 1 << a
                handled.update(((u, v), (v, u)))
            else:
                parent = row_image.get(u, u)
                for _ in range((length - 1) // 2):
                    x, y = state(v), state()
                    rows[parent] |= 1 << x
                    rows[x] |= 1 << y
                    rows[y] |= 1 << x
                    parent = y
                rows[parent] |= 1 << v
                handled.add((u, v))
    return tuple(rows), tuple(cycles), tuple(neighborhoods), tuple(projection), owner


def complete_basis(instance, selected):
    """Replace selected loop fillers in basis source cycles by the full basis.

    Matching robustness is preserved because all removed elements are
    matroid loops and all additions are old basis elements. Source-cycle
    quotas imply that the replacement does not increase cardinality.
    """
    selected = set(selected)
    if not selected <= set(range(len(instance.rows))):
        raise ValueError('Invalid state index.')
    original_size = len(selected)
    for j, cycle in enumerate(instance.source_cycles[:-1]):
        group = set(cycle)
        if len(selected & group) < 2:
            raise ValueError('A basis source-cycle quota is not satisfied.')
        pair = {j, instance.jobs + j}
        assert all(instance.neighborhoods[v] == 0 for v in group - pair)
        selected.difference_update(group - pair)
        selected.update(pair)
    assert instance.basis <= selected and len(selected) <= original_size
    return frozenset(selected)


def lift_covers(instance, first, second=None):
    """Lift covers of the two copies; reuse the first cover by default."""
    covers = tuple(map(set, (first, first if second is None else second)))
    full = (1 << instance.jobs) - 1
    for chosen in covers:
        union = 0
        for i in chosen:
            union |= instance.sets[i]
        if union != full:
            raise ValueError('Both chosen option families must cover the universe.')
    return instance.basis | {instance.options[c][i] for c in (0, 1) for i in covers[c]}


def extract_covers(instance, selected):
    """Recover two covers with total size at most |selected|-2*|U|.

    The selected placement must be robust. The helper additionally checks
    source-cycle quotas and resulting coverage, not full controllability.
    """
    selected = set(selected)
    original_size = len(selected)
    selected = complete_basis(instance, selected)
    covers = (set(), set())
    for v in selected - instance.basis:
        if v in instance.owner:
            c, i = instance.owner[v]
            covers[c].add(i)
    for chosen in covers:
        union = 0
        for i in chosen:
            union |= instance.sets[i]
        if union != (1 << instance.jobs) - 1:
            raise ValueError('The placement cannot repair every completed-basis failure.')
    assert sum(map(len, covers)) <= original_size - len(instance.basis)
    return tuple(map(frozenset, covers))
