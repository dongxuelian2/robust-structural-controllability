"""Exact mixed-failure spanning for direct sums of truncated unit-partition
matroids, using the independence oracle alone.

PROMISE: each connected component (after discarding loops) is a rank
truncation of a unit-capacity partition matroid. The connectivity
decomposition is recovered from fundamental circuits of a basis.
"""
from truncated_partition_fast import solve as solve_component


def solve(weights, independent, k, vulnerable=None):
    weights = tuple(weights)
    n = len(weights)
    if type(k) is not int or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if any(type(w) is not int or w < 0 for w in weights):
        raise ValueError("weights must be nonnegative integers")
    ground = frozenset(range(n))
    vulnerable = ground if vulnerable is None else frozenset(vulnerable)
    if not vulnerable <= ground:
        raise ValueError("vulnerable elements must belong to the ground set")

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        a, b = find(x), find(y)
        if a != b:
            parent[b] = a

    basis = []
    for e in range(n):
        if independent(frozenset(basis + [e])):
            basis.append(e)
    basis_set = frozenset(basis)

    for e in range(n):
        if e in basis_set:
            continue
        if not independent(frozenset({e})):
            continue  # A loop is its own rank-zero component.
        for b in basis:
            # b is in the fundamental circuit of e relative to the basis
            # if and only if exchanging b for e gives another basis.
            if independent((basis_set - {b}) | {e}):
                union(e, b)

    groups = {}
    for e in range(n):
        groups.setdefault(find(e), []).append(e)

    optimum = 0
    selected = set()
    infos = []
    for comp in groups.values():
        ids = tuple(comp)
        local_vulnerable = frozenset(i for i, x in enumerate(ids)
                                      if x in vulnerable)

        def local_independent(s):
            return independent(frozenset(ids[i] for i in s))

        value, witness, metadata = solve_component(
            [weights[x] for x in ids], local_independent, k, local_vulnerable
        )
        if value is None:
            return None, None, {"components": len(groups),
                                "infeasible": ids}
        optimum += value
        selected.update(ids[i] for i in witness)
        infos.append((len(ids), metadata["rank"], value))

    return optimum, frozenset(selected), {"components": len(groups),
                                         "ranks_and_costs": infos}
