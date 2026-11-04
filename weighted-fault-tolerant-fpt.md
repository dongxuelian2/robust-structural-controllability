# Weighted fault-tolerant spanning is FPT in rank and failure budget

## Question and result

Bentert, Fomin, Golovach and Morelle, [*Fault-Tolerant Matroid Bases*,
ESA 2025, Section 7](https://doi.org/10.4230/LIPIcs.ESA.2025.83), explicitly
ask whether their rank-plus-failure-budget FPT algorithm extends to weights
and to safe/vulnerable elements. The proof below addresses both questions.
It extends their uniform-set and bounded-solution arguments; the essential
change is weight-ordered selection and retaining lower-rank branches **even
when the uniform-set size threshold is reached**.

The source's [full version](https://arxiv.org/html/2506.22010v1) still lists
these questions. Targeted public searches on 4 November 2026 found no later
resolution. This identifies a specific published open question, not a
guarantee that no subsequent result exists.

**Theorem.** Given a rank-r matroid M on n elements by an independence
oracle, nonnegative binary-encoded integer weights, a vulnerable subset V,
and k>=0, a minimum-weight set B satisfying

    rank(B minus F)=r for every F subset B intersect V, |F|<=k

can be found, or infeasibility certified, in f(r,k) times a polynomial in
the input size and oracle-call cost. Setting V=E gives weighted ordinary
fault-tolerant spanning. The algorithm constructs an FPT-sized subset
containing an optimum; it is not a polynomial-size or standard encoded
kernel claim. It does not improve the general actuator-placement
approximation factor, and does not impose source-SCC quotas.

Safe elements are never deleted. The displayed selected-set condition is
the intended nonuniform model, rather than testing the full ground set
after failure. Loops may be discarded. If r=0 the empty set is optimal.
Assume r>=1 below. Write

    q=(k+1)r, H=binom(q-1,r-1), t_h=(h-1)H+q.

## A small optimum exists, also with weights and safe elements

For any feasible B, successively take a maximal independent subset I_j
of the remaining elements, for j=0,...,k, and let B0 be their union.
Its size is at most q. For any allowed failure F subset B0 of size at most k,
one of these k+1 disjoint layers avoids F. Every omitted element belongs
to the closure of that layer, by maximality at the time it was selected.
Consequently B minus F is contained in cl(B0 minus F). Feasibility of B
therefore implies feasibility of B0. Nonnegative weights give w(B0)<=w(B).
Starting with a minimum-weight solution thus yields an optimum of size<=q.
This argument needs only a bound on the number of failures, not that every
element is vulnerable.

## Weight-respecting exchange lemma

Let A be a subset containing elements of only one failure status, with
rank(A)=h>=1. Suppose Z subset A is h-uniform: every h-element subset
of Z is independent, and rank(Z)=h. Suppose |Z|>=t_h. Let B be any
feasible set with |B|<=q, and let x be in (B intersect A) minus Z with
w(z)<=w(x) for all z in Z.

Then there is y in Z minus B such that B'=(B minus {x}) union {y} is
feasible, has the same size, and has weight no larger than B.

Proof. Put Y=B minus {x}. For every allowed F subset Y intersect V of
size<=k for which rank(Y minus F)<r, feasibility of B implies that this
rank is exactly r-1 and x is outside the hyperplane cl(Y minus F).
Every such hyperplane has a basis of r-1 elements from Y; hence there
are at most binom(q-1,r-1)=H distinct hyperplanes. Each contains at most
h-1 elements of Z. Otherwise h independent elements of Z would span A,
forcing x into the hyperplane. At most q-1 elements of Z belong to B.
Since |Z|>=(h-1)H+q, some y in Z minus B avoids all these hyperplanes.
If there are no deficient failure sets, the same selection is immediate.

For an allowed failure F of B', if y is not failed, then F subset Y and
y repairs any deficient rank. If y is failed, it is vulnerable, so x is
also vulnerable. Replace y in F by x to obtain a legal failure of B;
the surviving set in both cases is Y minus (F minus {y}), which spans.
Thus B' is feasible. The weight and cardinality conclusions follow. QED.

## Ordered recursive core

The procedure Core(A) below preserves the following strong property:
for every feasible B with |B|<=q, there is a feasible B' of the same size
and no larger weight such that

    B' minus A = B minus A, and B' intersect A subset Core(A).

Here A has a single failure status and contains no loops. It need not be
a flat. Let h=rank(A).

1. If A is empty, return empty. If h=1, return the q cheapest elements
   of A, or all of A if fewer. If |A|<=t_h, returning A is also valid.
2. Otherwise scan A in nondecreasing weight order, breaking ties by a
   fixed element order. Begin with Z empty. Accept x exactly when every
   subset of Z union {x} of size at most h is independent. Stop upon
   accepting t_h elements or exhausting A.
3. For every (h-1)-subset T of Z, recursively compute Core(A_T), where
   A_T=A intersect cl(T). Return Z union all these child cores.

Step 3 is performed in **both** stopping cases. Tests in step 2 only
need h-subsets containing x once |Z|>=h-1; before that a single
independence test suffices. The final Z has at least h elements: before
reaching h its rule is ordinary greedy independence and A has rank h.
Thus it is h-uniform and spans A. Every child has rank exactly h-1.

**Preservation proof.** Induct on h. For h=1, all elements of A are
nonloop parallel elements. Any x outside the q cheapest can be exchanged
for an unused one of those q elements with the same status and no larger
weight. Parallel exchange preserves ranks and allowed failures. Repeat.

For h>1 let U be the union of the child sets A_T. Every rejected scanned
element x belongs to U: rejection exhibits a dependent set containing x
and at most h-1 retained elements; since the retained elements are
independent up to size h, x is spanned by them. Pad that set to an
(h-1)-subset of final Z. Also Z subset U.

Therefore, if the scan was exhausted, A=U. If the size threshold was met,
every x in A minus U was unscanned, so all of Z is no more expensive than
x. Apply the exchange lemma successively to elements of B in A minus U.
Every replacement lies in Z subset U and removes one element outside U.
The resulting feasible set has its A-part in U, the same size and no
greater weight. Its outside-A part is unchanged.

Now process each child A_T, applying its induction hypothesis to the
current solution. Each operation changes only elements in that child
and introduces elements only in that child's core. All child cores are
subsets of the final output. Thus an element not in the final output,
once removed, can never be reintroduced. Since every initial element of
the A-part lies in at least one child, after all children are processed
the A-part is contained in the output. This proves the strong property.

Apply the procedure separately to V and E minus V after removing loops.
The two ground subsets are disjoint, so the two compressions do not
interfere. Their union W contains an optimum of size at most q.

## Bound and exact algorithm

Define K_0=0, K_1=q and, for 2<=h<=r,

    K_h=t_h+binom(t_h,h-1)K_(h-1).

Every rank-h call returns at most K_h elements. This includes the
small-input shortcut, because t_h<=K_h. The sequence is nondecreasing,
so |W|<=2K_r (or K_r if all elements are vulnerable).
The number of recursive nodes, retained elements per node, and subsets
tested at a node are functions of r,k. Each node scans at most n elements
and constructs its children by polynomially many rank/independence
queries times a parameter-only factor. Ranks can themselves be computed
greedily from an independence oracle. Sorting and binary-weight
comparisons have polynomial bit cost. Thus core construction is FPT.

Enumerate all subsets B of W with |B|<=q. For each, enumerate allowed
failures of size at most k and test the surviving rank. Choose the
minimum-weight feasible set. Both enumerations have parameter-only size.
If none is feasible, the optimum-preservation property certifies
infeasibility. Nonnegative rational weights work as well with exact
arithmetic; the reference code takes nonnegative integers.

## Direct application to strongly connected state graphs

For a strongly connected state digraph, let T be the transversal matroid
of right-state vertices matchable by state edges and M=T*. Its rank is
the matching deficiency d. With one dedicated actuator per selected
state, matching feasibility after an actuator failure is precisely
spanning in M after deleting that selected element. The state vertices
and state edges themselves do not fail. This uses the standard
[dedicated-input matching and accessibility criterion, Theorem 3 of
Pequito, Kar and Aguiar](https://pages.up.pt/~up519521/publications/TAC2016_StrSys.pdf).

If d>0, a spanning surviving set is nonempty, so strong connectivity
supplies accessibility automatically. The theorem therefore gives an
exact weighted algorithm parameterized by d+k for this actuator model,
including safe and vulnerable actuators. If d=0, only accessibility
remains: choose the cheaper of one cheapest safe actuator and k+1
cheapest vulnerable actuators, omitting unavailable options. General
multi-source-SCC accessibility quotas require an additional argument
and are not covered by this corollary.

## Why simply sorting the old stopping rule is insufficient

In rank two with k=1, take two parallel classes of size two and weights
0,0 and 1,1, followed by five distinct additional directions of weight
100 each. This is a rank-two truncation of a partition matroid. Here q=4
and t_2=7. A greedy 2-uniform set reaches the threshold but retains only
one element from each cheap class. Within that set a robust solution
costs at least 101; the original optimum uses both cheap pairs and costs
2. The lower-rank branches retain them. The important difference is not
sorting alone: the recursion must survive the successful threshold case.

A higher-rank version uses a rank-three truncation: three parallel pairs
with costs (0,1) each, followed by 79 singleton classes of costs 13,...,91.
The optimum is 3, from all three cheap pairs; a uniform-only output has
optimum 13. Here t_3=26. The implemented recursion retains 29 of the 85
elements, including all cheap pairs. For rank-three truncations under one
failure, the exact optimum is the minimum of the four cheapest distinct
class minima and the three cheapest class-pair sums, so this test needs
no enumeration of all 85-element subsets.

## Verification

`weighted_fault_core.py` implements the core and the exact enumeration
using an arbitrary rank oracle. Run the two independent audits with
Python 3.10+ and no third-party packages:

    python audit_weighted_core.py
    python audit_weighted_independent.py

The first compares against every selected subset and every allowed
failure in its small-instance suites, including exhaustive binary
weights on specified rank-two, graphic, sparse-paving and nonrepresentable
Vamos matroids, all status assignments on a specified nine-element
instance, and seeded wider-weight tests. It also exercises successful
threshold termination together with rejected cheap elements.
The recorded run passed 18,779 comparisons, including 3,751 cases with
both threshold termination and rejection, and 42 infeasible cases.
The second uses independent binary-linear rank tables and repeated
worst-case-deletion dynamic programming on 1,200 seeded instances with
1,...,12 elements, plus the 85-element exact threshold example above.
Finite checks supplement the theorem; no linear representation is
assumed in the proof.
