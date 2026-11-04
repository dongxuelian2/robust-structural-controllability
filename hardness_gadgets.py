"""Cost-preserving subcubic reduction from Vertex Cover; Python 3.10+."""
from dataclasses import dataclass


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


def vertex_cover_graph(edges):
    """Return a connected, loopless, bipartite, bidirected subcubic instance.

    Source vertices are nonnegative integers; isolated vertices can be discarded.
    For m>0 source edges, target OPT = 2*m + minimum vertex cover size.
    No actuator is charged for a port or a component connector.
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
