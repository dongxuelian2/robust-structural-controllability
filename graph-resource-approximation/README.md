# Robust graph-resource assignment: path covers and approximation

This directory develops the simple graph-resource subclass of unit-cost, one-actuator-failure structural controllability, including the full four-source-state hardness construction in [PR #8](https://github.com/dongxuelian2/robust-structural-controllability/pull/8). It does **not** improve the general arbitrary-state-graph 3/2 guarantee in [PR #4](https://github.com/dongxuelian2/robust-structural-controllability/pull/4) or establish its optimality.

## Model and results

There is one job per vertex of a simple undirected graph. Each ordinary resource serves its two endpoints; each singleton resource serves one vertex. Select resources so that every job can still be matched after deletion of any one selected resource. Costs are one. Singleton copies may be distinct resources at the same vertex; parallel ordinary resources are excluded from the 4/3 theorem.

- **Implemented 4/3 algorithm:** permits arbitrary available singleton resources and an arbitrary mandatory subset of them. Forbidden singletons are removed before solving. The solver uses canonical singleton choices, minimum degree cover via blossom matching, and cycle-component augmentation.
- **Exact path-cover identities:** when the ordinary graph has no isolated vertices, with singleton completion, optimum is `n + pc(G)`; with two mandatory singletons at distinct `s,t`, it is `n + pc_{s,t}(G)`. A related identity holds without singletons under the stated minimum-degree and component assumptions.
- **Theoretical 8/7 algorithm transfer:** for a graph without isolated vertices and with all completion singletons available, the Berman-Karpinski maximum-path-edge theorem and a two-pendant-vertex reduction give 8/7 with two mandatory singleton resources. The full four-source-state control construction also has optimum `n + pc_{s,t}(G)`, including all auxiliary actuator choices, so receives the same bound.

The 8/7 claim imports the cited algorithmic theorem. **The Berman-Karpinski approximation algorithm is not implemented here.** Exact enumerators validate the reductions on small instances; they are not polynomial implementations of 8/7.

## Proofs

- [4/3 proof and control construction](four-thirds-graph-resources.md): canonicalization, matching gadget, control interpretation, and implementation scope.
- [Path-cover identities and 8/7 transfer](bicircular-pathcover-equivalence.md): constructive maps, mandatory anchors, pendant normalization, full-control optimum, and primary source.

Primary imported theorem: Berman and Karpinski, *8/7-Approximation Algorithm for (1,2)-TSP (Extended Version)*, ECCC TR05-069 revision 2 (June 2006), Theorem 2 and Section 5. [Primary PDF](https://eccc.weizmann.ac.il/report/2005/069/revision/2/download). The required objective is maximum **path-edge count**, with guarantee at least 6/7 of optimum, not minimum path-cover cardinality with a multiplicative 7/6 guarantee.

## Use the implemented solver

Python 3.10 or later; no third-party packages.

```python
from robust_assignment_solver import solve_graph_resources

n = 4
edges = [(0, 1), (1, 2), (2, 3)]
loop_vertices = [0, 1, 2, 3]
selected = solve_graph_resources(n, edges, loop_vertices, F=0)
resources = edges + [(v,) for v in loop_vertices]
print([resources[i] for i in range(len(resources)) if selected >> i & 1])
```

Resource bit positions list ordinary edges first, then singleton entries in the supplied order. `F` is a bit mask containing only mandatory singleton positions. The solver returns the selected-resource bit mask, or `None` if infeasible. Forbid a resource by removing it from the input and updating indices. Input validation rejects parallel ordinary edges. `polynomial_degree_cover.py` contains the matching gadget and blossom routine.

## Reproduce every check

From the repository root, enter the publication subdirectory and run:

```sh
cd graph-resource-approximation
python verify_four_thirds.py
python verify_control_four_thirds.py
python audit_bicircular_pathcover.py
python audit_anchored_pathcover.py
python audit_uncompleted_pathcover.py
python audit_pendant_anchors.py
python audit_control_pathcover_independent.py
```

All modules have guarded entry points; importing helpers does not run enumerations.

Recorded successful coverage:

- 1,896 mandatory-loop fixtures; 94,264 selected-resource subsets checked by independent matching. Includes exhaustive optimum, canonical-loop and polynomial degree-cover comparisons.
- 18 actual control graphs; 198 source-pattern outputs; 3,365 direct state/input matching and accessibility checks; 26,203 robust supports in exact-optimum enumerations.
- 901 singleton-completed, 499 two-anchor, and 423 uncompleted path-cover instances.
- Separate assignment-matching audit: 259 full source-copy/quota fixtures and 34,416 candidate supports, independently checking the full-control optimum and the path-cover lift.
- 259 pendant-anchor instances; 10,220 augmented path covers; 703 covers satisfying the cited 6/7 edge guarantee, all passing the resulting 8/7 inequality.

The checks are finite supplements to the proofs. No publication-priority, current-best, general inapproximability, or approximation-optimality claim is made.
