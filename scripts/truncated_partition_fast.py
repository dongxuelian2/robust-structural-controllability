"""Faster polynomial exact solver for mixed safe/vulnerable failures on
rank truncations of unit-capacity partition matroids (independence oracle).

Uses SMAWK on a Monge extension of the convex min-plus transition.
The accompanying proof note states the matroid-class promise.
"""


def _leftmost_row_minima(size, value):
    """SMAWK on a square totally monotone matrix accessed through value."""
    minima = [0] * size

    def search(rows, columns):
        if not rows:
            return
        reduced = []
        for column in columns:
            while reduced:
                row = rows[len(reduced) - 1]
                if value(row, reduced[-1]) <= value(row, column):
                    break
                reduced.pop()
            if len(reduced) < len(rows):
                reduced.append(column)

        search(rows[1::2], reduced)
        offset = 0
        for index in range(0, len(rows), 2):
            row = rows[index]
            right = minima[rows[index + 1]] if index + 1 < len(rows) else reduced[-1]
            best_column = reduced[offset]
            best_value = value(row, best_column)
            while reduced[offset] < right:
                offset += 1
                column = reduced[offset]
                current = value(row, column)
                if current < best_value:
                    best_column, best_value = column, current
            minima[row] = best_column

    search(range(size), range(size))
    return minima


def _vulnerable_transition(dp, prefix, inf):
    maximum = len(prefix) - 1

    def value(z, u):
        count = z - u
        if count < 0:
            cost = -inf * count
        elif count > maximum:
            cost = prefix[-1] + inf * (count - maximum)
        else:
            cost = prefix[count]
        return dp[u] + cost

    minima = _leftmost_row_minima(len(dp), value)
    return [(min(value(z, u), inf), u) for z, u in enumerate(minima)]


def solve(weights, independent, k, vulnerable=None):
    weights = tuple(weights)
    if type(k) is not int or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if any(type(w) is not int or w < 0 for w in weights):
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
        return 0, frozenset(), {"rank": 0, "threshold": 0}
    k_eff = min(k, len(vulnerable))
    inf = sum(weights) + 1
    options = []
    for cls in classes:
        ss = sorted((x for x in cls if x not in vulnerable),
                    key=lambda x: (weights[x], x))
        vv = sorted((x for x in cls if x in vulnerable),
                    key=lambda x: (weights[x], x))
        options.append((ss[0] if ss else None, vv))

    best_cost = best_set = best_tau = None
    for tau in range(1, k_eff + 2):
        Q = (r - 1) * tau + k_eff + 1
        dp = [inf] * (Q + 1)
        dp[0] = 0
        histories = []
        for safe, vuln in options:
            m = min(tau, len(vuln))
            prefix = [0]
            for x in vuln[:m]:
                prefix.append(prefix[-1] + weights[x])

            new_dp = [inf] * (Q + 1)
            back = [None] * (Q + 1)
            for z, (cost, predecessor) in enumerate(_vulnerable_transition(dp, prefix, inf)):
                if cost < inf:
                    new_dp[z] = cost
                    back[z] = (predecessor, False, z - predecessor)
                if safe is not None:
                    u = max(0, z - tau)
                    val = min(dp[u] + weights[safe], inf)
                    if val < new_dp[z]:
                        new_dp[z] = val
                        back[z] = (u, True, 0)
            dp = new_dp
            histories.append(back)

        cost = dp[Q]
        if cost >= inf or (best_cost is not None and cost >= best_cost):
            continue
        chosen = set()
        z = Q
        for i in range(len(classes) - 1, -1, -1):
            old, use_safe, count = histories[i][z]
            safe, vuln = options[i]
            if use_safe:
                chosen.add(safe)
            else:
                chosen.update(vuln[:count])
            z = old
        assert z == 0
        best_cost, best_set, best_tau = cost, frozenset(chosen), tau

    return best_cost, best_set, {"rank": r, "classes": len(classes),
                                 "threshold": best_tau,
                                 "effective_failures": k_eff}
