# Weighted fault tolerant matroid spanning

- [Read the paper (PDF)](weighted-fault-tolerant-spanning.pdf).
- [Edit the LaTeX source](weighted-fault-tolerant-spanning.tex).
- [Read the shorter proof note](../weighted-fault-tolerant-fpt.md).

The manuscript proves exact FPT for minimum-weight selected-set spanning with parameter rank plus failure budget, including safe and vulnerable elements. It preserves the explicit distinction from the literal nonuniform formulation in ESA 2025. It is a manuscript for review, not a submitted publication; the author list is unset.

## Compile

From the repository root, with a standard TeX Live or MiKTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/weighted-fault-tolerant-spanning.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/weighted-fault-tolerant-spanning.tex
```

## Verify

The existing root-level files are the accompanying implementation and checks; no duplicate copies are needed:

- [weighted_fault_core.py](../weighted_fault_core.py): rank-oracle core construction and exact solver.
- [audit_weighted_core.py](../audit_weighted_core.py): exhaustive subset and failure checks on the specified finite suites.
- [audit_weighted_independent.py](../audit_weighted_independent.py): independent deletion-DP checks and targeted counterexamples.

With Python 3.10 or newer and no third-party packages, run from the repository root:

```sh
python audit_weighted_core.py
python audit_weighted_independent.py
```

The reference solver implements the general core-and-enumeration algorithm, not the optional polynomial-time shortcuts for rank one and zero failures described in the paper. Finite checks supplement the proofs.

The 5 November refinement uses the classical Bollobas set-pairs inequality to sharpen the deficient-hyperplane bound to `binom(r+k-1,k)`. The 85-element example now retains 15 elements rather than 29, with the same exact optimum. The independent audit includes direct checks of the counting lemma, separate from the solver.

