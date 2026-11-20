# Weighted resilient spanning with arbitrary disjoint cardinality quotas

An exact algorithm for weighted resilient spanning with disjoint purchase quotas. The companion [robust-accessibility theorem](mixed-accessibility-quotas.md) treats group survival under mixed statuses.

## Problem and claim

Input: rank-r finite matroid M=(E,I) via independence oracle; nonnegative binary-encoded weights w; safe/vulnerable partition; failure budget k; pairwise-disjoint groups S_1,...,S_s with integer cardinality lower bounds beta_i >= 1. Seek minimum w(B) subject to rank(B minus F)=r for every F subset of selected vulnerable elements of size at most k, and |B intersect S_i| >= beta_i.

Groups with beta_i>|S_i| imply infeasibility. Zero quotas can be discarded. The groups may contain loops. The cardinality constraints apply to the purchased set. In the all-vulnerable case robust source-SCC accessibility is a cardinality quota beta_i=k+1.

Claim: for r,k>=1, this problem is solvable in

    2^{O(b log b + br log H)} poly(n,L),
    b=(k+1)r, H=binom(r+k-1,k),

independent of the number and sizes of quota groups in the exponential factor. The r=0 case is separate independent cheapest-quota completion. At k=0 the same construction gives FPT with b=r and H=1.

## Quota completion cost

Let G0 be the elements outside the groups. For any forced set C define g(C) as the minimum weight of a superset satisfying all quotas. It is computed by including C, then adding the cheapest max(0,beta_i-|C intersect S_i|) elements of S_i minus C, independently for each group. Outside the groups nothing else is added. Nonnegative weights justify this choice. Fix deterministic weight tie-breaking.

For each group let a_i be the beta_i-th smallest original weight in S_i, including loops and both statuses. Define adjusted weight

    v(e)=max(0,w(e)-a_i) if e in S_i; v(e)=w(e) if e in G0.

Lemma A (marginal completion).
For x not in Y, g(Y+x)-g(Y)>=v(x). If Y intersect S_i is empty and x in S_i, equality holds. If x in G0 equality always holds.

Proof: In a group where Y already has at least beta_i elements, the marginal is w(x)>=v(x). Otherwise put t=beta_i-|Y intersect S_i|>0 and let c be the t-th cheapest weight in S_i minus Y. We have c<=a_i: among the beta_i original cheapest elements at most |Y intersect S_i| were removed, leaving at least t. The marginal of forcing x is max(0,w(x)-c). This formula is valid whether x is among the current t cheapest or not, including ties. Since c<=a_i, it is at least max(0,w(x)-a_i). If Y avoids the group then c=a_i.

Lemma B (same-group exchange).
If x,y are distinct, neither in Y, and both belong to the same group or both G0, with w(y)<=w(x), then g(Y+y)<=g(Y+x).

Proof for a group: choose an optimum completion D containing Y+x. If y is already in D, D is a feasible completion of Y+y. If not, D-x+y retains the same group counts and has no larger weight. For G0 the claim is additive. Note D need not be resilient; this lemma concerns completion cost only.

Ordering: process elements by lexicographic (v(e),w(e),fixed ID). Within any one group this order refines nondecreasing original weight. This remains true for tied v=0. The same is true for G0.

## Small resilient witness

Every resilient B contains a resilient subset C of size at most b, by the standard k+1 disjoint maximal-independent-layer argument, also with safe/vulnerable statuses. Removing loops is allowed in C but not in the quota ground set.

Consequently the desired optimum equals min g(C) over resilient loopless C with |C|<=b. In one direction, an optimal full B has such a C and g(C)<=w(B); in the other, any quota completion of resilient C remains resilient by monotonicity of survivor spanning and has cost g(C).

## Flat-saturation core with color caps

Apply the following procedure separately to safe and vulnerable NONLOOP elements, always using the ambient matroid's flats. Each group S_i is a color, and G0 is an additional color. Set L0=A0=0 and for h>=1 define

    Lh = H L(h-1)+b,
    Ah = H A(h-1)+b Lh+b.

Process elements in the ordering above. Maintain W and all constraints

    |W intersect Q intersect color i| <= L(rank Q),
    |W intersect Q| <= A(rank Q)

for every ambient flat Q and every color i. Reject x if a containing flat is saturated in either the TOTAL constraint or x's OWN COLOR constraint; otherwise insert it. Only x's own color count can change, so this rule preserves all invariants. Each safety-status core has size <=Ar.

Implementation: every saturated flat Q of rank h is spanned by its retained elements (respectively its retained own-color elements). Otherwise their closure P has rank d<h and the SAME relevant retained count, contradicting strict increase of Lh/Ah and the invariant for P. Thus enumerate independent subsets I of W of size at most r; Q=cl(I). For each candidate x count W intersect Q and W intersect Q intersect color(x), and test x in Q. This detects all witnesses. Rank-zero flats contain no scanned elements. Each scan has parameter-only work times n; a direct bound O(n(r+1)(Ar+1)^{r+1}) oracle queries suffices per safety status.

## Exchange proof

Fix any resilient C with |C|<=b and x in C outside the union core. Let Y=C-x. The standard deficient-hyperplane argument gives at most H distinct hyperplanes P_j=cl(Y-F), for allowed vulnerable failures F of Y that leave rank r-1, and all P_j avoid x. A replacement y of the same failure status, outside C and outside every P_j, makes Y+y resilient. This is exactly the exchange predicate of the competing preprint and the existing repository proof; if y fails map that failure to x.

Take the saturated witness flat Q at the time x was rejected, rank h. Every retained witness element precedes x in adjusted order and is of the same failure status. Since Q contains x while P_j does not, Q intersect P_j is a flat of rank at most h-1.

Case 1, own-color saturation: Z=W intersect Q intersect color(x) has Lh elements. Each P_j contains at most L(h-1) of them, so at least b avoid all P_j. At most b-1 are in Y. Choose another y. Within the same color earlier ordering implies w(y)<=w(x); Lemma B yields g(Y+y)<=g(Y+x).

Case 2, total saturation: Z=W intersect Q has Ah elements. Exclude G0 and every group containing an element of Y. There are at most |Y|+1<=b forbidden colors, costing at most b Lh elements. The union of P_j costs at most H A(h-1) further elements. At least

    Ah - b Lh - H A(h-1) = b

remain. Every remaining y belongs to a group untouched by Y, so y not in Y. It is not x because x was rejected. It has v(y)<=v(x). Lemma A gives

    g(Y+y)-g(Y)=v(y)<=v(x)<=g(Y+x)-g(Y).

Thus this replacement is resilient and cannot increase completion cost. Both witness types reduce the number of selected elements outside the final union core by one. Repeating places an optimum small witness in that core.

## Size and search complexity

Writing sums with H^0=1,

    Lr = b sum_{j=0}^{r-1} H^j,
    Ar = b sum_{j=0}^{r-1} H^j + b^2 sum_{j=0}^{r-1}(j+1) H^j.

(The convolution follows from Ah=H A(h-1)+b Lh+b.) Hence Ar <= b r H^{r-1}+b^2 r(r+1)H^{r-1}/2 for H>=1. For r,k>=1, H>=2 except r=1, and log(2Ar+1)=O(log b+r log H); for r=1 Ar=b+b^2 so this also holds.

Enumerate subsets C of union core of size <=b. Test resilience by all <=2^b selected-vulnerable failure patterns and greedy rank queries. Compute g(C) by sorting group weights once and completing quotas, polynomial in input size. Choose cheapest completed set. If no feasible witness exists, report infeasibility.

The search has exponent O(b log b+br log H). Core construction's exponent is O((r+1)(log b+r log H)); r+1<=b for k>=1 absorbs it in the same bound. The number of colors only contributes polynomial bookkeeping/input cost, not a parameter exponent. Exact sums and comparisons retain polynomial bit complexity.

## Relation to flat saturation

Kundu, Nandi and Saurabh give the same exponent for the unquotaed problem. Their Lemma 5 supplies the single-cap saturation rule. The quota extension uses adjusted completion prices and separate per-color and total caps. [*When Elements Fail: Weighted Fault-Tolerant Matroid Bases*, Theorem 3.1 and Lemma 5](https://arxiv.org/html/2611.16333v1).

The companion [robust-accessibility theorem](mixed-accessibility-quotas.md) charges an additional parameter for groups with a cheaper safe option and a finite vulnerable backup. It gives an algorithm on arbitrary state digraphs parameterized by matching deficiency, failure budget, and the number of such groups.

## Reproduction

Run `python scripts/audit_weighted_quotas.py` with Python 3 and assertions enabled. The deterministic run checks 96,797 marginals, 71,782 same-group inequalities, 179 optimum comparisons, and 874 sampled exchanges. It exercises 6 color-cap and 772 total-cap discards. The small-instance tests include loops, zero weights, and mixed statuses.
