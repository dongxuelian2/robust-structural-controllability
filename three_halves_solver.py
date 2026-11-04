"""Unit-cost single-actuator-failure 3/2 approximation, no third-party dependencies.

Uses the original exact quota relaxation and a maximum general-graph matching.
Graph-operation complexity O(n^3), storage O(n^2); bitmasks carry normal Python
big-integer costs. maximum_matching is an Edmonds augmenting-forest algorithm
with blossom contraction. No greedy maximal-matching substitute is valid.
"""
from collections import deque
from two_approx_solver import solve as quota_solver


def maximum_matching(adjacency):
    """Return symmetric mate array for a maximum matching of a simple graph.

    Input is a sequence of symmetric neighbor collections without self loops.
    This function does not mutate its input. Each augmenting search contracts
    odd alternating cycles using base representatives, then lifts the path via
    its updated predecessor links. O(v^3) time and O(v^2) input storage.
    """
    n = len(adjacency)
    graph = [tuple(sorted(set(ns))) for ns in adjacency]
    if any(not isinstance(v, int) or v < 0 or v >= n or v == u
           for u, ns in enumerate(graph) for v in ns):
        raise ValueError('Expected simple undirected graph adjacency.')
    neighbor_sets = [set(ns) for ns in graph]
    if any(u not in neighbor_sets[v] for u, ns in enumerate(graph) for v in ns):
        raise ValueError('Expected symmetric adjacency.')
    mate = [-1] * n

    def augment(root):
        base = list(range(n))
        parent = [-1] * n
        outer = [False] * n
        queue = deque([root])
        outer[root] = True

        def ancestor(u, v):
            seen = [False] * n
            while True:
                u = base[u]
                seen[u] = True
                if mate[u] == -1:
                    break
                u = parent[mate[u]]
            while True:
                v = base[v]
                if seen[v]:
                    return v
                v = parent[mate[v]]

        def mark_arm(v, common, child, contracted):
            while base[v] != common:
                contracted[base[v]] = True
                contracted[base[mate[v]]] = True
                parent[v] = child
                child = mate[v]
                v = parent[mate[v]]

        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if base[u] == base[v] or mate[u] == v:
                    continue
                if v == root or (mate[v] != -1 and parent[mate[v]] != -1):
                    common = ancestor(u, v)
                    contracted = [False] * n
                    mark_arm(u, common, v, contracted)
                    mark_arm(v, common, u, contracted)
                    for x in range(n):
                        if contracted[base[x]]:
                            base[x] = common
                            if not outer[x]:
                                outer[x] = True
                                queue.append(x)
                elif parent[v] == -1:
                    parent[v] = u
                    if mate[v] == -1:
                        x = v
                        while x != -1:
                            previous = parent[x]
                            following = mate[previous] if previous != -1 else -1
                            mate[x] = previous
                            if previous != -1:
                                mate[previous] = x
                            x = following
                        return True
                    partner = mate[v]
                    outer[partner] = True
                    queue.append(partner)
        return False

    for root in range(n):
        if mate[root] == -1:
            augment(root)
    return tuple(mate)


def solve(rows):
    """Return robust placement and inspectable matching/repair certificates."""
    rows = tuple(rows)
    initial = quota_solver(rows)
    if initial is None:
        return None
    n = len(rows)
    P = initial['relaxation_actuators']
    D = initial['critical_states']
    base = sum(1 << b for b in P)
    incoming = [[] for _ in range(n)]
    for u, row in enumerate(rows):
        for v in range(n):
            if row >> v & 1:
                incoming[v].append(u)
    left_mate = [-1] * n
    for u, v in initial['matching']:
        left_mate[u] = v
    covers = {e: [] for e in range(n) if not base >> e & 1}
    for index, b in enumerate(D):
        reached_right = {b}
        reached_left = set()
        queue = deque([b])
        while queue:
            v = queue.popleft()
            for u in incoming[v]:
                if u in reached_left:
                    continue
                reached_left.add(u)
                e = left_mate[u]
                assert e != -1, 'A fragile point cannot reach a free left vertex.'
                if e not in reached_right:
                    reached_right.add(e)
                    queue.append(e)
        for e in reached_right - {b}:
            covers[e].append(index)
    adjacency = [set() for _ in D]
    witness = {}
    individual = {}
    for e, covered in covers.items():
        for i in covered:
            individual.setdefault(i, e)
            for j in covered:
                if i < j:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
                    witness.setdefault((i, j), e)
    mates = maximum_matching(adjacency)
    repairs = set()
    pair_certificates = []
    singleton_certificates = []
    for i, j in enumerate(mates):
        if j == -1:
            e = individual[i]
            repairs.add(e)
            singleton_certificates.append((D[i], e))
        elif i < j:
            e = witness[(i, j)]
            repairs.add(e)
            pair_certificates.append((D[i], D[j], e))
    chosen = sorted(set(P) | repairs)
    q = len(D) - len(pair_certificates)
    return {
        'actuators': tuple(chosen),
        'relaxation_actuators': P,
        'relaxation_value': len(P),
        'critical_states': D,
        'repair_states': tuple(sorted(repairs)),
        'matching': initial['matching'],
        'repair_supports': tuple((e, tuple(D[i] for i in covered))
                                 for e, covered in covers.items()),
        'pair_repairs': tuple(pair_certificates),
        'singleton_repairs': tuple(singleton_certificates),
        'pair_matching_size': len(pair_certificates),
        'repair_budget': q,
        'opt_lower_bound': max(len(P), 2 * q),
    }
