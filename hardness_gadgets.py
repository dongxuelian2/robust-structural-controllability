"""Cost-preserving subcubic reduction from Vertex Cover; Python 3.10+."""
from dataclasses import dataclass
from collections import Counter


@dataclass
class Instance:
    rows: tuple
    labels: tuple
    edges: tuple
    centers: frozenset
    forced_pairs: tuple
    repair_sets: tuple
    terminals: dict
    owner: dict
    # Subdivision vertices project to parallel copies of these original states.
    projection: tuple


def vertex_cover_graph(edges, port_order=None):
    """Return a connected, loopless, bipartite, bidirected subcubic instance.

    Source vertices are nonnegative integers; isolated vertices can be discarded.
    For m>0 source edges, target OPT = 2*m + minimum vertex cover size.
    No actuator is charged for a port or a component connector.
    For planarity, supply port_order[v] as the cyclic neighbor order of a
    planar source embedding. The default ordering need not preserve planarity.
    """
    edges = tuple(sorted({tuple(sorted(e)) for e in edges}))
    if not edges or any(len(e) != 2 or e[0] == e[1] or
                        any(not isinstance(v, int) or v < 0 for v in e)
                        for e in edges):
        raise ValueError('Use a nonempty simple graph on nonnegative integers.')
    incidence = {}
    for j, edge in enumerate(edges):
        for v in edge:
            incidence.setdefault(v, []).append(j)
    if port_order is not None:
        edge_index = {edge: j for j, edge in enumerate(edges)}
        for v, incident in incidence.items():
            ordered = tuple(port_order[v])
            expected = {w for j in incident for w in edges[j] if w != v}
            if len(ordered) != len(expected) or set(ordered) != expected:
                raise ValueError('Each port order must list exactly the neighbors.')
            incidence[v] = [edge_index[tuple(sorted((v, w)))] for w in ordered]
    labels, rows, centers = [], [], set()

    def state(label, center=False):
        index = len(rows)
        labels.append(label)
        rows.append(0)
        if center:
            centers.add(index)
        return index

    def link(a, b):
        rows[a] |= 1 << b
        rows[b] |= 1 << a

    pairs, elements = [], []
    for j, _ in enumerate(edges):
        c = state(('c', j), True)
        a, b, x = (state((name, j)) for name in ('a', 'b', 'x'))
        for t in (a, b, x):
            link(c, t)
        pairs.append((a, b))
        elements.append(x)
    q, ports, terminals, owner = {}, {}, {}, {}
    for v, incident in sorted(incidence.items()):
        previous = None
        q[v] = []
        for i, j in enumerate(incident):
            p = state(('p', v, i), True)
            t = state(('q', v, i))
            link(p, t)
            link(p, elements[j])
            if previous is not None:
                link(previous, p)
            ports[v, j] = p
            q[v].append(t)
            owner[t] = v
            previous = t
        terminals[v] = q[v][-1]

    # A prefix pair provides a degree-one center in each source component.
    # Joining these centers changes only freely matchable center columns.
    unseen = set(incidence)
    connectors = []
    while unseen:
        root = min(unseen)
        todo, component = [root], {root}
        for v in todo:
            for j in incidence[v]:
                for w in edges[j]:
                    if w not in component:
                        component.add(w)
                        todo.append(w)
        unseen.difference_update(component)
        p0 = state(('connector', root), True)
        q0 = state(('prefix', root))
        link(p0, q0)
        link(q0, ports[root, incidence[root][0]])
        connectors.append(p0)
    for a, b in zip(connectors, connectors[1:]):
        link(a, b)

    repair = []
    for j, edge in enumerate(edges):
        candidates = {elements[j]}
        owner[elements[j]] = edge[0]
        for v in edge:
            i = incidence[v].index(j)
            candidates.update(q[v][i:])
        repair.append(frozenset(candidates))
    return Instance(tuple(rows), tuple(labels), edges, frozenset(centers),
                    tuple(pairs), tuple(repair), terminals, owner,
                    tuple(range(len(rows))))


def lift_cover(instance, cover):
    return frozenset(v for pair in instance.forced_pairs for v in pair) | {
        instance.terminals[v] for v in cover}


def extract_cover(instance, actuators):
    """For a feasible unsubdivided placement, return a cover of size <= |S|-2m.

    This helper rejects failed necessary constraints. Robustness is independently
    checked by the verifier, not assumed on the strength of this helper.
    """
    actuators = set(actuators)
    forced = {v for pair in instance.forced_pairs for v in pair}
    if not forced <= actuators or any(not (actuators & r)
                                      for r in instance.repair_sets):
        raise ValueError('The placement violates a necessary robust constraint.')
    return frozenset(instance.owner[v] for v in actuators - forced
                     if v in instance.owner)


def protect_coloops(rank, ground, support):
    """Complete a spanning support using <= one addition per initial coloop.

    `rank` accepts sets; the ambient matroid must be coloop-free. For a
    robust multiset of parallel copies, this produces a distinct original
    robust set of no larger cost. Accessibility is a separate condition.
    """
    ground, support = set(ground), set(support)
    d = rank(ground)
    if rank(support) != d:
        raise ValueError('The projected support does not span.')
    critical = [v for v in support if rank(support - {v}) < d]
    for v in critical:
        removed = support - {v}
        if rank(removed) == d:
            continue
        replacement = next((w for w in sorted(ground - support)
                            if rank(removed | {w}) == d), None)
        if replacement is None:
            raise ValueError('The ambient matroid has a coloop.')
        support.add(replacement)
    return frozenset(support)


def normalize_subdivision(base, subdivided, actuators):
    """Map a robust subdivided placement to a no-larger original placement.

    `base` is the unsubdivided instance. `subdivided` is obtained by any
    number of calls to subdivide_ports(base). Feasibility is a precondition;
    the audit checks the actual actuator-failure model independently.
    """
    forced = {v for pair in base.forced_pairs for v in pair}
    relevant = forced | set(base.owner)
    counts = Counter(subdivided.projection[v] for v in set(actuators)
                     if subdivided.projection[v] in relevant)
    label_index = {label: v for v, label in enumerate(base.labels)}
    for j, (a, b) in enumerate(base.forced_pairs):
        ca, cb = counts[a], counts[b]
        if ca + cb < 2:
            raise ValueError('Fewer than two selected copies of a forced pair.')
        if ca >= 2 and cb >= 2:
            x = label_index['x', j]
            counts[x] = max(1, counts[x])
        counts[a] = counts[b] = 1
    result = frozenset(v for v, count in counts.items() if count)
    assert len(result) <= len(set(actuators))
    return result


def subdivide_ports(instance):
    """Replace each center--target edge by a length-three bidirected path.

    Center--center connector bridges are retained. Every new target is a
    parallel copy, in the dual matching matroid, of its incident old target;
    new centers are matroid loops. Degrees and planarity are preserved, and
    lengths of all cycles are multiplied by three.
    """
    rows = [0] * len(instance.rows)
    labels = list(instance.labels)
    centers = set(instance.centers)
    projection = list(instance.projection)

    def link(a, b):
        rows[a] |= 1 << b
        rows[b] |= 1 << a

    for u, row in enumerate(instance.rows):
        for v in range(u + 1, len(instance.rows)):
            if not row >> v & 1:
                continue
            if u in instance.centers and v in instance.centers:
                link(u, v)
                continue
            c, t = (u, v) if u in instance.centers else (v, u)
            a, b = len(rows), len(rows) + 1
            rows.extend((0, 0))
            labels.extend((('copy', t, a), ('subcenter', c, b)))
            projection.extend((instance.projection[t], -1))
            centers.add(b)
            link(c, a)
            link(a, b)
            link(b, t)
    return Instance(tuple(rows), tuple(labels), instance.edges,
                    frozenset(centers), instance.forced_pairs,
                    instance.repair_sets, instance.terminals, instance.owner,
                    tuple(projection))
