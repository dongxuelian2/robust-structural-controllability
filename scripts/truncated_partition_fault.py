"""Exact weighted safe/vulnerable fault-tolerant spanning on rank-r
truncations of unit-capacity partition matroids (independence-oracle input).

PRECONDITION: after discarding loops, a set is independent exactly when
it has at most one element per nonloop parallel class and size at most r.
The class promise is not recognized or verified on arbitrary oracle inputs.
Weights must be nonnegative integers; the budget k may be arbitrarily large.
"""


def solve(weights, independent, k, vulnerable=None):
    """Return (optimum cost, selected IDs, certificate), or (None, None, ...)."""
    weights = tuple(weights)
    if type(k) is not int or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if any(type(c) is not int or c < 0 for c in weights):
        raise ValueError("weights must be nonnegative integers")
    n = len(weights)
    ground = frozenset(range(n))
    vulnerable = ground if vulnerable is None else frozenset(vulnerable)
    if not vulnerable <= ground:
        raise ValueError("vulnerable elements must belong to the ground set")

    classes = []
    for x in range(n):
        if not independent(frozenset({x})):
            continue
        for cls in classes:
            if not independent(frozenset({cls[0], x})):
                cls.append(x)
                break
        else:
            classes.append([x])

    basis = []
    for cls in classes:
        x = cls[0]
        if independent(frozenset(basis + [x])):
            basis.append(x)
    r = len(basis)
    if r == 0:
        return 0, frozenset(), {"rank": 0, "classes": 0, "threshold": 0}
    if len(classes) < r:
        raise AssertionError("rank cannot exceed number of classes")

    # A failure budget greater than the total vulnerable ground size is
    # equivalent to the effective budget k_eff <= n.
    k_eff = min(k, len(vulnerable))
    class_options = []
    for cls in classes:
        ss = sorted((x for x in cls if x not in vulnerable),
                    key=lambda x: (weights[x], x))
        vv = sorted((x for x in cls if x in vulnerable),
                    key=lambda x: (weights[x], x))
        class_options.append((ss, vv))

    best_cost = best_set = best_tau = None
    for tau in range(1, k_eff + 2):
        target = (r - 1) * tau + k_eff + 1
        options_by_class = []
        for ss, vv in class_options:
            options = [(0, 0, ())]
            if ss:
                options.append((tau, weights[ss[0]], (ss[0],)))
            running = 0
            for j, x in enumerate(vv[:tau], start=1):
                running += weights[x]
                options.append((j, running, tuple(vv[:j])))
            options_by_class.append(options)

        dp = [None] * (target + 1)
        dp[0] = 0
        trace = []
        for options in options_by_class:
            next_dp = [None] * (target + 1)
            back = [None] * (target + 1)
            for old_cov, old_cost in enumerate(dp):
                if old_cost is None:
                    continue
                for choice, (cov, cost, _) in enumerate(options):
                    coverage = min(target, old_cov + cov)
                    val = old_cost + cost
                    if next_dp[coverage] is None or val < next_dp[coverage]:
                        next_dp[coverage] = val
                        back[coverage] = (old_cov, choice)
            dp = next_dp
            trace.append(back)
        cost = dp[target]
        if cost is None or (best_cost is not None and cost >= best_cost):
            continue
        selected = set()
        cov = target
        for i in range(len(classes) - 1, -1, -1):
            old_cov, choice = trace[i][cov]
            selected.update(options_by_class[i][choice][2])
            cov = old_cov
        assert cov == 0
        best_cost, best_set, best_tau = cost, frozenset(selected), tau

    return best_cost, best_set, {"rank": r, "classes": len(classes),
                                 "threshold": best_tau,
                                 "effective_failures": k_eff}
