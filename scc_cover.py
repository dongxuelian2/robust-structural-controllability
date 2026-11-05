"""Set Cover via source-SCC quotas, with OPT = 2*(|U| + tau)."""
from dataclasses import dataclass


@dataclass
class CoverInstance:
    rows: tuple
    sets: tuple
    jobs: int
    forced: frozenset
    options: tuple
    neighborhoods: tuple
    owner: dict


def cover_graph(n, sets, sparse=True):
    """Return a loopless bipartite control instance and its cover maps.

    Two copies of the covering instance give OPT=2*(n+tau), without parity
    assumptions or a special case for tau=1. In the sparse version every
    underlying undirected degree is two or three, except jobs with source
    frequency greater than two. Its maximum directed outdegree is three.
    """
    sets = tuple(sets)
    full = (1 << n) - 1
    if (n < 1 or not sets or any(s <= 0 or s & ~full for s in sets) or
            any(not any(s >> j & 1 for s in sets) for j in range(n))):
        raise ValueError('Use nonempty subsets covering a nonempty universe.')
    q = len(sets)
    count = 4 * n + 2 * q
    rows, neighborhoods, color = [0] * count, [0] * count, [0] * count
    forced, owner = frozenset(range(2 * n)), {}
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
    assert all(color[u] != color[v] for u, row in enumerate(rows)
               for v in range(len(rows)) if row >> v & 1)
    return CoverInstance(tuple(rows), sets, n, forced, options,
                         tuple(neighborhoods), owner)


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
    return instance.forced | {instance.options[c][i] for c in (0, 1) for i in covers[c]}


def extract_covers(instance, selected):
    """Recover two covers with total size at most |selected|-2*|U|.

    The selected placement must be robust. The helper additionally checks
    the forced states and the resulting coverage, not full controllability.
    """
    selected = set(selected)
    if not instance.forced <= selected or not selected <= set(range(len(instance.rows))):
        raise ValueError('The placement omits forced states or has invalid indices.')
    covers = (set(), set())
    for v in selected - instance.forced:
        if v in instance.owner:
            c, i = instance.owner[v]
            covers[c].add(i)
    for chosen in covers:
        union = 0
        for i in chosen:
            union |= instance.sets[i]
        if union != (1 << instance.jobs) - 1:
            raise ValueError('The placement cannot repair every forced-basis failure.')
    assert sum(map(len, covers)) <= len(selected) - len(instance.forced)
    return tuple(map(frozenset, covers))
