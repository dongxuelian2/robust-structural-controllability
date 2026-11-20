# Robust Structural Controllability under Actuator Failures

Research on dedicated actuator placement under failures, including
approximation algorithms, exact algorithms, and fault-tolerant matroid
spanning. The state graph defines a transversal matroid; its dual translates
matching feasibility into a spanning condition. Source strongly connected
components impose additional accessibility quotas.

## Results and models

- The [3/2 approximation](paper/robust-dedicated-actuator-placement.md)
  runs in cubic time on arbitrary directed state graphs for unit-cost
  dedicated actuators at distinct states under one actuator failure.
- [Factor-two approximation](two-approximation.md): for a nonempty directed
  state graph, unit-cost dedicated actuators at distinct states, and one
  actuator failure, a quota-first algorithm gives a factor-two guarantee and
  runs in cubic graph operations. The note also gives a family where the
  independent-repair strategy approaches factor two.
- [Bounded-core algorithm](bounded-core-theorem.md): the same one-failure
  placement model has an exact algorithm polynomial for each fixed matching
  deficiency `d`, obtained by enumerating robust cores of at most `2d` states.
- [Disjoint-quota theorem](quota-core-proof.md): unit-cost robust matroid
  spanning with disjoint quotas is fixed-parameter tractable in matroid rank
  and failure budget. Applied to dedicated actuators, the parameters are
  matching deficiency and the number of failures.
- [Weighted fault-tolerant spanning](paper/weighted-fault-tolerant-spanning.pdf)
  ([LaTeX source](paper/weighted-fault-tolerant-spanning.tex),
  [proof note](weighted-fault-tolerant-fpt.md)): weighted selected-set
  spanning in oracle matroids is fixed-parameter tractable in rank and
  failure budget, including the stated safe/vulnerable extension.

The dedicated placement results use one actuator per selected state and unit
costs. The weighted matroid theorem concerns selected-set spanning. The
[flat-saturation revision](https://github.com/dongxuelian2/robust-structural-controllability/pull/16)
gives the updated comparison with the quartic-exponent oracle time bound.

The starting reference for robust structural controllability is Diogo Poças,
Sérgio Pequito, and Guilherme Ramos, [*On the NP-hardness of the robust
minimal structural controllability problem*](https://doi.org/10.1016/j.automatica.2026.113074),
*Automatica* 190 (2026), 113074.

## Verification

Python 3.10 or newer and the standard library are sufficient. From the
repository root:

```sh
python scripts/verify_core_theorem.py
python scripts/check_two_approx.py
python scripts/audit_weighted_core.py
python scripts/audit_weighted_independent.py
```

The checks compare algorithms with exhaustive or independent finite-instance
oracles. Each proof note specifies its test domain and corresponding solver.
The [paper guide](paper/README.md) gives manuscript build instructions and
links to the integrated control manuscript.
