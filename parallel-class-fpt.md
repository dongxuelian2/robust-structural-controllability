# Fault-tolerant spanning with few parallel classes

## Result and scope

Let M be an n-element matroid given by an independence oracle, with nonnegative
integer costs, safe and vulnerable elements, and an integer failure budget k.
A feasible selected set must span M after every deletion of at most k of its
selected vulnerable elements. Let p be the number of nonloop parallel classes.

**Theorem.** A minimum-cost feasible set, or a report of infeasibility,
can be computed in f(p) times a polynomial in the input encoding length and n,
with polynomial degree independent of p and k. Oracle calls are counted at unit
cost. In particular, for every fixed finite field F_Q, this problem is FPT in
rank r alone on F_Q-representable matroids. The finite-field consequence needs
only an independence oracle and the representability promise, not a supplied
representation. If Q varies, the parameter is (r,Q), not r alone.

The proof combines the standard cocircuit characterization, parallel-class
compression, and fixed-dimensional integer programming. Kundu, Nandi and
Saurabh [1, Section 4.5] already use parallel classes to retain a cheapest safe
element and k+1 cheapest vulnerable elements per class. Their stated subsequent
search bound is exponential in k. Here class multiplicities are integer
variables, removing k from the parameter function. This is a restricted-model
parameter improvement, **not** a uniform improvement of their general oracle
bound or arbitrary-field algorithm. The parameter dependence obtained via
integer programming can be very large; for small k their algorithm may be
better. Publication priority for this combination has not been established.

## Hyperplane characterization

Write H for a hyperplane of M. A selected set B is feasible exactly when every H
has either a selected safe element outside H, or at least k+1 selected vulnerable
elements outside H. Indeed, if neither holds, deleting all selected elements
outside H is an allowed failure and leaves a nonspanning set. Conversely, any
nonspanning survivor set lies in a hyperplane, whose selected outside elements
must all have failed. This characterization is standard; see [1, Section 2].

Rank zero is immediate: the empty set is optimal. Otherwise discard loops.
Replace k by min(k, |V|), which changes no admissible effective failure.
For each parallel class C, let a_C be its cheapest safe cost, when available.
Sort its vulnerable costs as w_C,1 <= ... <= w_C,m. There is an optimum in which:

1. A class with a selected safe element has exactly one such element, of cost
   a_C, and no selected vulnerable element.
2. Other classes contain at most k+1 selected vulnerable elements, chosen as
   the cheapest elements of that class.

For (1), a surviving safe element spans every other element of its class, so
removing the others preserves every survivor rank. For (2), retaining k+1
copies guarantees that class survives any allowed failure; deleting excess
copies cannot affect feasibility. Exchanging identities within the same class
and safety type preserves feasibility by the characterization. Nonnegative
costs ensure these operations do not increase cost.

## An integer program with 3p variables

Put m_C = min(number of vulnerable elements in C, k+1) and
P_C(j) = sum_{i=1}^j w_C,i, with P_C(0)=0. Use integer variables s_C, x_C, t_C.
The variable s_C records the choice of one cheapest safe element; x_C counts
selected vulnerable elements; t_C represents their prefix cost. Impose

    0 <= s_C <= 1,
    0 <= x_C <= m_C,
    x_C <= (k+1)(1-s_C),
    t_C >= 0.

If C has no safe element, impose s_C=0 and use a_C=0 in the objective.
For every hyperplane H impose

    sum_{C outside H} [x_C + (k+1)s_C] >= k+1.                 (1)

Parallel classes are wholly inside or outside a flat. Since s_C is binary,
(1) is exactly the hyperplane characterization, including when no failures
occur. It is not a probabilistic or fractional relaxation.

For every 0 <= j < m_C impose

    t_C >= P_C(j) + w_C,j+1 (x_C-j).                          (2)

Minimize sum_C (a_C s_C + t_C). The sorted slopes make the piecewise-linear
interpolation of P_C convex. Each line in (2) supports that interpolation;
for an integer x_C in [0,m_C], at least one line attains P_C(x_C), except
m_C=0 when t_C>=0 suffices. Hence the smallest permitted t_C is exactly
P_C(x_C). All these values are integers. This proves exact equality of the
integer-program optimum and the original optimum. Recover the selected set
using its prescribed cheapest elements.

## Oracle construction and bit complexity

Discover loops and parallel classes using singleton and pair independence
queries, in O(n^2) queries. Pick one representative per class and obtain rank
by greedy independence tests. Enumerate subsets of these p representatives.
For every subset of rank r-1, compute its closure on the representatives and
keep the distinct closures. This lists every hyperplane of the simplification:
each has a basis among the representatives. It takes 2^p poly(p) oracle calls
and produces at most 2^p constraints (1).

There are 3p integer variables and O(n+p+2^p) constraints. Prefix sums and
intercepts in (2) have O(log(n+1)+log(W_max+1)) bits; coefficients involving k
have O(log(k+1)) bits. Thus the integer-program encoding length is bounded by
a function of p times a polynomial in the original encoding length.

Fixed-dimensional integer programming [2] solves feasibility in f(d)L^{O(1)}
bit operations for d integer variables, with an absolute polynomial degree.
For optimization, add an integer cost-budget inequality and binary-search the
integer interval from zero to the sum of all input costs. The number of calls
is polynomial in the input bit length. A feasible solution with all t_C at
their prefix values has cost at most that sum, so the initial budget also
detects infeasibility. We may bound each t_C by this sum without losing such a
solution. This proves the claimed parameterized bound, including recovery of
an optimum. The proof invokes the fixed-dimensional algorithm, not a claim
that the audit script implements it.

For rank r over F_Q there are at most (Q^r-1)/(Q-1) one-dimensional subspaces,
hence at most that many parallel classes. Substituting this bound proves the
corollary. The corresponding statement for arbitrary fields does not follow:
there may be unboundedly many directions at fixed rank. In fact [3, Theorem 4]
proves hardness already at fixed rank at least three over the rationals.

## Reproducible finite checks

Run `python audit_parallel_class_fpt.py` with Python 3.10 or newer. No packages
are needed. It compares exhaustive subset-and-failure optimization against
independent class-count enumeration and ambient-hyperplane inequalities over
GF(2) and GF(3). It also checks the epigraph identity at every integer prefix
length. Seeded cases are reproducible, not exhaustive over all matroids;
within each fixture all selected subsets and allowed failures are enumerated.
The test is a check of the reduction, not an implementation or empirical speed
claim for Lenstra's algorithm, and not a substitute for the proof.

## References

[1] M. Kundu, S. Nandi, S. Saurabh, *When Elements Fail: Weighted Fault-Tolerant
Matroid Bases*, arXiv:2611.16333v1, 12 November 2026, Sections 2 and 4.5.
https://arxiv.org/html/2611.16333v1

[2] H. W. Lenstra, Jr., *Integer Programming with a Fixed Number of Variables*,
Mathematics of Operations Research 8(4), 538-548, 1983.
https://doi.org/10.1287/moor.8.4.538
Author copy: https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1983i/art.pdf

[3] M. Bentert, F. V. Fomin, P. A. Golovach, L. Morelle, *Fault-Tolerant Matroid
Bases*, ESA 2025, Theorem 4. This theorem also already gives a polynomial-time
weighted algorithm at rank at most two in the all-vulnerable model.
https://doi.org/10.4230/LIPIcs.ESA.2025.83
