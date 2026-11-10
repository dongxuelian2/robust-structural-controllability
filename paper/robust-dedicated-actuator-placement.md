# Robust dedicated actuator placement manuscript

**Preprint manuscript:** The LaTeX source and checked-in PDF reflect the current editorial revision. The preprint has not been submitted to a journal.

[PDF](robust-dedicated-actuator-placement.pdf) · [LaTeX](robust-dedicated-actuator-placement.tex)

The manuscript gives complete proofs of three results:

- A cubic-time 3/2 approximation for unit-cost, distinct-state dedicated actuators under one actuator failure on arbitrary directed graphs.
- Weighted selected-set fault-tolerant spanning in an independence-oracle matroid, FPT in rank plus failure budget, with safe and vulnerable elements and the set-pairs hyperplane bound.
- Exact unit-cost spanning with arbitrarily many disjoint quota-(k+1) groups, FPT in rank plus failure budget. The control corollary applies to arbitrary digraphs independently of the source-SCC count.

The quota algorithm completes its small robust core using original elements, including loops. It is not a bounded full-solution kernel. Weighted quotas and the optimal polynomial-time approximation factor remain open.

## Compile

Use a standard TeX Live or MiKTeX installation:

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
```

## Reproduce the algorithms

All implementations use Python 3.10 or newer and the standard library.

The [weighted and quota snapshot](https://github.com/dongxuelian2/robust-structural-controllability/tree/9d62c26633643d05973e3ecefc5d36e7a8b426c1) contains `quota_fault_core.py`, `weighted_fault_core.py`, their audits, and the shared `bounded_core_solver.py`. These files are unchanged by the manuscript revision.

The 3/2 implementation from PR #4 is fixed at [commit b510683](https://github.com/dongxuelian2/robust-structural-controllability/tree/b510683764600ccd318eb7db566c35a816c2edc2). That snapshot contains the [solver](https://github.com/dongxuelian2/robust-structural-controllability/blob/b510683764600ccd318eb7db566c35a816c2edc2/three_halves_solver.py), [proof note](https://github.com/dongxuelian2/robust-structural-controllability/blob/b510683764600ccd318eb7db566c35a816c2edc2/three-halves-approximation.md), both audits below, and its `two_approx_solver.py` and `bounded_core_solver.py` dependencies. No branch-tip version or merge of PR #4 is required.

For isolated checkouts of the two fixed versions:

```sh
git clone https://github.com/dongxuelian2/robust-structural-controllability.git
cd robust-structural-controllability
git worktree add --detach ../robust-quota-check 9d62c26633643d05973e3ecefc5d36e7a8b426c1
git worktree add --detach ../robust-approx-check b510683764600ccd318eb7db566c35a816c2edc2
(cd ../robust-quota-check && python audit_quota_core.py && python audit_weighted_core.py && python audit_weighted_independent.py)
(cd ../robust-approx-check && python audit_blossom.py && python audit_three_halves_independent.py)
```

The audits compare the implementations with independent exact calculations on small matroids and directed graphs, and include larger threshold and exhaustion fixtures. Mathematical proofs are in the manuscript; finite tests check implementation behavior.

Sparse-hardness research remains separate in [PR #3](https://github.com/dongxuelian2/robust-structural-controllability/pull/3).
