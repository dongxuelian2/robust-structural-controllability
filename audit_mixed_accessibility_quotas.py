"""Finite checks of prices, completion, and robust group compression."""
from math import inf
from random import Random

from audit_weighted_quotas import core, gf2rank, resilient, subsets


def group_options(weights, colors, safe):
    return [
        (
            sorted(weights[e] for e in range(len(weights)) if colors[e] == group and safe[e]),
            sorted(weights[e] for e in range(len(weights)) if colors[e] == group and not safe[e]),
        )
        for group in range(1, max(colors, default=0) + 1)
    ]


def group_prices(weights, colors, safe, required):
    prices = list(weights)
    protected = set()
    for group, (secure, exposed) in enumerate(group_options(weights, colors, safe), 1):
        safe_price = secure[0] if secure else inf
        backup_price = sum(exposed[:required]) if len(exposed) >= required else inf
        if safe_price == backup_price == inf:
            return None, None
        if safe_price < backup_price < inf:
            protected.add(group)
            continue
        for element, color in enumerate(colors):
            if color != group:
                continue
            if safe[element]:
                prices[element] -= min(safe_price, backup_price)
            elif backup_price < inf:
                prices[element] = max(0, weights[element] - exposed[required - 1])
    return prices, protected


def completion(forced, weights, colors, safe, required):
    cost = sum(weights[e] for e in forced)
    for group in range(1, max(colors, default=0) + 1):
        if any(colors[e] == group and safe[e] for e in forced):
            continue
        missing = required - sum(colors[e] == group and not safe[e] for e in forced)
        if missing <= 0:
            continue
        safe_price = min((weights[e] for e in range(len(weights))
                          if colors[e] == group and safe[e]), default=inf)
        available = sorted(weights[e] for e in range(len(weights))
                           if colors[e] == group and not safe[e] and e not in forced)
        backup_price = sum(available[:missing]) if len(available) >= missing else inf
        cost += min(safe_price, backup_price)
    return cost


def groups_survive(selected, colors, safe, k):
    exposed = [e for e in selected if not safe[e]]
    groups = set(colors) - {0}
    return all(groups <= {colors[e] for e in set(selected) - set(failed)}
               for failed in subsets(exposed, k))


def best_witness(ground, weights, colors, safe, rank, r, k):
    return min((completion(set(chosen), weights, colors, safe, k + 1)
                for chosen in subsets(sorted(ground), (k + 1) * r)
                if resilient(chosen, rank, r, k, safe)), default=inf)


def main():
    generator = Random(20261120)
    marginal_checks = optimum_checks = protected_groups = 0
    for trial in range(160):
        n = generator.randrange(1, 10)
        weights = [generator.randrange(9) for _ in range(n)]
        colors = [generator.randrange(min(n, 3) + 1) for _ in range(n)]
        renumber = {g: i for i, g in enumerate(sorted(set(colors) | {0}))}
        colors = [renumber[g] for g in colors]
        safe = [generator.randrange(2) == 0 for _ in range(n)]
        vectors = [generator.randrange(8) for _ in range(n)]

        def rank(elements):
            return gf2rank(vectors[e] for e in elements)

        r = rank(range(n))
        k = generator.randrange(1, 3)
        prices, protected = group_prices(weights, colors, safe, k + 1)
        exact = min((sum(weights[e] for e in candidate)
                     for candidate in subsets(list(range(n)), n)
                     if groups_survive(candidate, colors, safe, k)
                     and resilient(candidate, rank, r, k, safe)), default=inf)
        optimum_checks += 1
        if prices is None:
            assert exact == inf
            continue
        protected_groups += len(protected)
        if r:
            retained, _ = core(rank, r, k, weights, colors, safe, {},
                               prices=prices, protected=protected)
        else:
            retained = set()
        assert best_witness(retained, weights, colors, safe, rank, r, k) == exact
        for base in subsets(list(range(n)), n):
            base = set(base)
            baseline = completion(base, weights, colors, safe, k + 1)
            for e in set(range(n)) - base:
                if colors[e] in protected:
                    continue
                delta = completion(base | {e}, weights, colors, safe, k + 1) - baseline
                assert delta >= prices[e]
                if colors[e] == 0 or all(colors[i] != colors[e] for i in base):
                    assert delta == prices[e]
                marginal_checks += 1

    assert protected_groups > 0
    assert [completion(chosen, [5, 4, 4], [1, 1, 1], [True, False, False], 2)
            for chosen in (set(), {1}, {2}, {1, 2})] == [5, 8, 8, 8]

    colors = [i // 4 + 1 for i in range(48)]
    safe = [i % 4 == 0 for i in range(48)]
    weights = [1 if i < 16 else 7 if safe[i] else 2 for i in range(48)]

    def rank_one(elements):
        return int(bool(elements))

    prices, protected = group_prices(weights, colors, safe, 2)
    retained, discards = core(rank_one, 1, 1, weights, colors, safe, {},
                              prices=prices, protected=protected)
    assert protected == {1, 2, 3, 4} and all(discards.values())
    assert best_witness(retained, weights, colors, safe, rank_one, 1, 1) == best_witness(
        set(range(48)), weights, colors, safe, rank_one, 1, 1
    )
    print(f'PASS: {optimum_checks} full optimum comparisons and {marginal_checks} price checks.')
    print(f'PASS: {protected_groups} protected groups, declining marginal, and both discard rules.')


if __name__ == '__main__':
    main()
