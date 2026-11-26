# Algorithms for Robust Dedicated Actuator Placement

[PDF](robust-dedicated-actuator-placement.pdf) · [LaTeX](robust-dedicated-actuator-placement.tex)

The revised arXiv manuscript centers on the matching–accessibility decomposition and proves three control results:

- A cubic-time 3/2 approximation for one actuator failure and unit costs on arbitrary directed graphs.
- A deterministic weighted algorithm with all locations vulnerable, fixed-parameter tractable in matching deficiency and failure budget using group-aware flat saturation.
- A one-sided randomized weighted algorithm for mixed safe and vulnerable locations on arbitrary digraphs using represented matroids and group-state coefficient compression.

The complete earlier standalone matroid results remain in the repository, including [weighted flat saturation](../weighted-flat-core.md), [weighted disjoint quotas](../weighted-disjoint-quotas.md), [mixed accessibility](../mixed-accessibility-quotas.md), and [represented matroids and sampling](../linear-source-quotas.md).

## Build

From the repository root, with TeX Live or MiKTeX:

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
pdflatex -interaction=nonstopmode -halt-on-error robust-dedicated-actuator-placement.tex
```

The manuscript is a single self-contained LaTeX source, including the figures and bibliography.

## Reproduce

From the repository root, with Python 3.10+ and the standard library:

```sh
python scripts/audit_three_halves_independent.py
python scripts/audit_weighted_quotas.py
python scripts/audit_linear_accessibility_quotas.py
```

The algorithms and proof notes are linked above. The sampled representation implementation uses an optimization-specific sampling bound; the paper gives the simpler full-representation union bound. The finite tests are implementation checks, not empirical industrial performance benchmarks.
