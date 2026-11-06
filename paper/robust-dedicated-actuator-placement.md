# Robust dedicated actuator placement manuscript

[PDF](robust-dedicated-actuator-placement.pdf) · [LaTeX](robust-dedicated-actuator-placement.tex)

The manuscript gives complete proofs of three results:

- A cubic-time 3/2 approximation for unit-cost, distinct-state dedicated actuators under one actuator failure on arbitrary directed graphs.
- Weighted selected-set fault-tolerant spanning in an independence-oracle matroid, FPT in rank plus failure budget, with safe and vulnerable elements and the sharp set-pairs hyperplane count.
- Exact unit-cost spanning with arbitrarily many disjoint quota-(k+1) groups, FPT in rank plus failure budget. The control corollary applies to arbitrary directed graphs without parameterizing the source-SCC count.

The quota core retains an optimum small robust core; completing its quotas can use original elements outside the core, including loops. This is not a bounded full-solution kernel or a weighted-quota theorem. The optimal approximation factor and publication priority remain unresolved. The manuscript includes AI-assistance disclosure.

Compile with an ordinary TeX Live installation:

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
```

From the repository root, using Python 3.10 or newer:

```sh
python audit_quota_core.py
python audit_weighted_core.py
python audit_weighted_independent.py
```

The new quota audit passes 23,020 matroid/quota comparisons and 1,750 independently computed control comparisons. It also checks a 43-state graph with 21 source SCCs and three threshold/exhaustion fixtures. The implementation reuses the existing `bounded_core_solver.py` for state matching and SCCs.

The 3/2 implementation and its blossom/full-rank audits remain in [PR #4](https://github.com/dongxuelian2/robust-structural-controllability/pull/4). They are not duplicated or implicitly merged here. That branch's `audit_blossom.py` and `audit_three_halves_independent.py` were separately rerun: 35,068 matching cases and 66,066 exhaustive plus 800 random state graphs passed. Finite audits validate implementation behavior, not theorem truth or priority.

The long sparse-hardness reduction chain remains separate in [PR #3](https://github.com/dongxuelian2/robust-structural-controllability/pull/3); no new numerical hardness bound is asserted in this manuscript.
