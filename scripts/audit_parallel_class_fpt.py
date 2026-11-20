"""Exact finite checks of parallel-class FPT reduction, not an ILP solver.

Python 3.10+, standard library only. Enumerate every subset and failure within
each small fixture, and independently enumerate canonical class choices.
"""
import itertools
import json
import random


def rank(vectors, prime):
    rows = [list(v) for v in vectors]
    if not rows:
        return 0
    pivot = 0
    for col in range(len(rows[0])):
        row = next((i for i in range(pivot, len(rows)) if rows[i][col]), None)
        if row is None:
            continue
        rows[pivot], rows[row] = rows[row], rows[pivot]
        inverse = pow(rows[pivot][col], -1, prime)
        rows[pivot] = [v * inverse % prime for v in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot:
                factor = rows[i][col]
                rows[i] = [(x - factor * y) % prime
                           for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def direction(vector, prime):
    first = next((v for v in vector if v), None)
    if first is None:
        return None
    return tuple(v * pow(first, -1, prime) % prime for v in vector)


def brute_optimum(vectors, safe, weights, k, prime):
    target = rank(vectors, prime)
    best = None
    for mask in range(1 << len(vectors)):
        chosen = [i for i in range(len(vectors)) if mask >> i & 1]
        vulnerable = [i for i in chosen if not safe[i]]
        feasible = True
        for count in range(min(k, len(vulnerable)) + 1):
            for failure in itertools.combinations(vulnerable, count):
                survivors = [vectors[i] for i in chosen if i not in failure]
                if rank(survivors, prime) != target:
                    feasible = False
                    break
            if not feasible:
                break
        if feasible:
            cost = sum(weights[i] for i in chosen)
            best = cost if best is None else min(best, cost)
    return best


def count_optimum(vectors, safe, weights, k, prime):
    if rank(vectors, prime) == 0:
        return 0
    # Test fixtures use a full-rank ambient representation in this branch.
    dimension = len(vectors[0])
    if rank(vectors, prime) != dimension:
        raise ValueError("Nonzero fixtures must have full ambient row rank")
    classes = {}
    for i, vector in enumerate(vectors):
        key = direction(vector, prime)
        if key is not None:
            classes.setdefault(key, []).append(i)
    keys = list(classes)
    options = []
    for key in keys:
        indices = classes[key]
        costs = sorted(weights[i] for i in indices if not safe[i])[:k + 1]
        prefix = [0]
        for cost in costs:
            prefix.append(prefix[-1] + cost)
        for x in range(len(prefix)):
            epigraph = max([0] + [prefix[j] + costs[j] * (x - j)
                                  for j in range(len(costs))])
            if epigraph != prefix[x]:
                raise AssertionError(("epigraph", costs, x))
        choices = [(0, x, cost) for x, cost in enumerate(prefix)]
        safe_costs = [weights[i] for i in indices if safe[i]]
        if safe_costs:
            choices.append((1, 0, min(safe_costs)))
        options.append(choices)
    # Every nonzero linear functional gives an ambient hyperplane. Duplicate
    # kernels are harmless. This construction does not use the brute solver.
    normals = [v for v in itertools.product(range(prime), repeat=dimension)
               if any(v)]
    outside = [[j for j, key in enumerate(keys)
                if sum(a * b for a, b in zip(key, normal)) % prime]
               for normal in normals]
    best = None
    for chosen in itertools.product(*options):
        if all(sum(chosen[j][1] + (k + 1) * chosen[j][0] for j in indices)
               >= k + 1 for indices in outside):
            cost = sum(entry[2] for entry in chosen)
            best = cost if best is None else min(best, cost)
    return best


def check(vectors, safe, weights, k, prime):
    expected = brute_optimum(vectors, safe, weights, k, prime)
    actual = count_optimum(vectors, safe, weights, k, prime)
    if actual != expected:
        raise AssertionError((prime, vectors, safe, weights, k, actual, expected))


def main():
    seed = 371
    rng = random.Random(seed)
    random_cases = 0
    for prime in (2, 3):
        for dimension in (1, 2, 3):
            for _ in range(50):
                n = rng.randrange(dimension, 10)
                vectors = [tuple(rng.randrange(prime) for _ in range(dimension))
                           for _ in range(n)]
                if rank(vectors, prime) != dimension:
                    continue
                safe = [rng.randrange(2) for _ in vectors]
                weights = [rng.randrange(7) for _ in vectors]
                k = rng.randrange(n + 3)
                check(vectors, safe, weights, k, prime)
                random_cases += 1
    fixtures = [
        ([], [], [], 0, 2),
        ([(0, 0), (0, 0)], [0, 1], [0, 7], 10, 2),
        ([(1,), (1,), (0,)], [0, 0, 1], [0, 3, 0], 2, 2),
        ([(1,), (1,), (1,)], [1, 0, 0], [8, 1, 2], 1, 2),
        ([(1, 0), (0, 1)], [1, 1], [0, 5], 10**20, 2),
        ([(1, 0), (0, 1)], [0, 0], [0, 5], 0, 2),
        ([(1, 0), (2, 0), (0, 1), (0, 2)],
         [1, 0, 0, 0], [8, 0, 1, 2], 1, 3),
    ]
    for fixture in fixtures:
        check(*fixture)
    print(json.dumps({"seed": seed, "random_full_rank_cases": random_cases,
                      "boundary_fixtures": len(fixtures),
                      "result": "all exact optima and epigraph values match"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
