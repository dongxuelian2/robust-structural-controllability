"""Faster polynomial exact solver for mixed safe/vulnerable failures on
rank truncations of unit-capacity partition matroids (independence oracle).

Uses a Monge min-plus convolution optimized by divide and conquer.
The matroid-class promise is required; see the accompanying proof note.
"""


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

            v_costs = [inf] * (Q + 1)
            v_arg = [0] * (Q + 1)

            def fill(lo, hi, opt_lo, opt_hi):
                if lo > hi:
                    return
                z = (lo + hi) // 2
                start = max(opt_lo, z - m)
                end = min(opt_hi, z)
                minval = None
                minu = None
                for u in range(start, end + 1):
                    val = dp[u] + prefix[z - u]
                    if minval is None or val < minval:
                        minval, minu = val, u
                assert minu is not None
                v_costs[z] = min(minval, inf)
                v_arg[z] = minu
                fill(lo, z - 1, opt_lo, minu)
                fill(z + 1, hi, minu, opt_hi)

            fill(0, Q, 0, Q)
            new_dp = [inf] * (Q + 1)
            back = [None] * (Q + 1)
            for z in range(Q + 1):
                new_dp[z] = v_costs[z]
                back[z] = (v_arg[z], False, z - v_arg[z])
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
