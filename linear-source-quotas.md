# Weighted dedicated actuator placement under mixed failures

Let \(M\) be a rank-\(r\) matroid represented by columns over a field.
Each element has a nonnegative binary-encoded cost and is designated
safe or vulnerable. For pairwise disjoint groups \(G_1,\ldots,G_g\),
the survivors of the purchased set must span \(M\) and meet every group
after any deletion of at most \(k\) purchased vulnerable elements.

**Represented-matroid theorem.** For \(r,k\ge1\), the exact optimum or
infeasibility can be found deterministically in
\[
 2^{O(kr(r+k)\log(r+1))}\operatorname{poly}(n,L)
\]
field operations and cost bookkeeping, independent of the number and
cost structure of the groups in the exponent. Here \(n\) is the number
of columns and \(L\) is the input bit length. Over the rationals the
same bound holds in bit operations. For \(r=0\) or \(k=0\) there are
polynomial-time algorithms.

**Actuator theorem.** Let a nonempty directed state graph allow at most
one dedicated actuator at each state, with nonnegative binary-encoded
actuator costs, and mixed safe/vulnerable failures. Write
\(d=n-\rho\), where \(\rho\) is the size of a maximum state matching.
For \(d,k\ge1\), a randomized algorithm runs in
\[
 2^{O(kd(d+k)\log(d+1))}
 \operatorname{poly}(n,L,\log(1/\delta))
\]
bit operations for error tolerance \(0<\delta<1\). Every returned
placement is feasible for the original graph. On a feasible input it
returns a true optimum with probability at least \(1-\delta\). On
an infeasible input every trial returns without a placement. The
parameter function is independent of the source-SCC count and of the
number \(u\) of exceptional groups in the
[oracle algorithm](mixed-accessibility-quotas.md).

The group counters and rank-checked sampling follow suggestions by
Saket Saurabh. The coefficient method is from Kundu, Nandi and
Saurabh [1, Sections 4.3–4.4].

## Small witnesses inside purchases

Use \(P\) for the whole purchase and \(C\subseteq P\) for its spanning
witness. If \(P\) survives the spanning requirement, choose \(k+1\)
disjoint layers by successively taking a maximal independent subset
of the remaining elements of \(P\). Their union \(C\) has size
at most \(m=(k+1)r\). For any allowed failure \(F\subseteq C\),
one layer avoids \(F\). Every element of \(P\setminus C\) is in the
closure of that layer, so \(C\setminus F\) spans \(P\setminus F\)
and hence spans \(M\). Thus every feasible purchase admits a small
resilient witness.

Conversely, a resilient witness guarantees the spanning condition
for all purchases containing it: for a failure of the whole purchase,
apply resilience to its intersection with \(C\). Group feasibility is
equivalent to purchasing a safe element or at least \(k+1\) vulnerable
elements in each group. Sufficiency follows from the failure budget,
and necessity follows by failing all of the group's purchases when
they are vulnerable and number at most \(k\). Loop columns can serve
as quota-only purchases. The cost charged is \(w(P)\).

## Coefficients for common future witnesses

Prescribe final safe and vulnerable witness counts with total at most
\(m\). Within a prefix state, let \(W\) be a partial witness of rank
\(p\), and fix its safe and vulnerable counts. There are \(q\) future
slots determined by these counts, with the safe slots first. A future
witness \(Z=(z_1,\ldots,z_q)\) assigns its elements of each status
to the corresponding slots in the fixed processing order.

Let \(T_j\) be a formal slot tag and \(X_h\) a formal \(r\)-vector
for \(0\le h<r-p+k\). For a fixed \(Z\), substitute
\[
 X_h\longmapsto \sum_{j=1}^q u_j T_j^h a_{z_j}
\]
in the polynomial ring, with independent variables \(u_j\) and with
the tags unchanged. Here \(a_e\) is the input column of element \(e\).
For a failed set \(J\) of vulnerable future slots, put
\[
 f_J(t)=\prod_{j\in J}(t-T_j)=\sum_\nu c_{J,\nu}t^\nu,
 \qquad
 Y_{J,h}=\sum_\nu c_{J,\nu}X_{h+\nu}.
\]
Under the substitution, \(Y_{J,h}\) becomes the moment of the
surviving columns weighted by \(u_jf_J(T_j)\).

For an allowed past failure \(I\subseteq W\cap V\), take a column
basis \(A_I\) of \(W\setminus I\), where \(V\) is the vulnerable
set. Put \(\eta=r-\operatorname{rank}(W\setminus I)\). For
\(|I|+|J|\le k\), use the factor
\[
 D_{I,J}=\det[\,A_I\mid Y_{J,0}\mid\cdots\mid Y_{J,\eta-1}\,],
\]
with value \(1\) for \(\eta=0\). Its largest moment index is
\(\eta-1+|J| \le r-p+k-1\). After substitution, the coefficient
of \(\prod_{j\in H}u_j\) for an \(\eta\)-subset \(H\) of surviving
future slots is the input minor for \(A_I\) and those slots,
multiplied by their Vandermonde determinant and
\(\prod_{j\in H} f_J(T_j)\). These last factors are nonzero.
Distinct choices \(H\) have distinct \(u\)-monomials. Thus the
factor evaluates to a nonzero polynomial exactly when the
survivors for this failure pair span \(M\). This applies over
every field, since the tags are distinct indeterminates. Over
\(\mathbb Q\) the tags may be fixed to distinct integers \(1,\ldots,q\).

Set \(\Phi_W=\prod_{|I|+|J|\le k}D_{I,J}\), where the pairs
range over the allowed past and future failures. The evaluation
ring is an integral domain, so \(\Phi_W\) evaluates to a nonzero
polynomial exactly when \(W\cup Z\) is a resilient witness. For
fixed \(Z\) the substitution is a single linear map on coefficient
vectors, common to the candidates in the state. This simultaneous
test supplies a common completion under all allowed failures.

## Dynamic program and complexity

Process the disjoint groups consecutively and the ungrouped elements
in a separate block. Skip an element, buy it for the witness, or buy
it for the group alone. Either purchase pays the full cost. The
state stores the partial witness counts and rank, and the current
group records a safe-purchase flag and the number of vulnerable
purchases capped at \(k+1\). At the group boundary require the flag
or a saturated count and reset the two counters. Check every
prescribed final pair of witness counts. Empty required groups
imply infeasibility.

Within a state, sort candidates by total purchase cost and retain
each coefficient vector outside the span of those already kept.
For any discarded candidate, its vector is in the span of retained
vectors of no greater cost. If its witness works with a fixed
future witness \(Z\), the common substitution is nonzero on at
least one of those cheaper vectors. The same future quota-only
purchases also work, since all completed groups passed and the
current group has identical counters. After a group reset, any
pooled candidates have passed that group.

For correctness, fix the future decisions of any feasible marked
purchase, including its witness labels. Inductively a retained
prefix of no greater cost admits those decisions. Taking the next
decision generates a successor before compression, in its actual
state. The compression argument preserves a successor for the
remaining decisions. At the end use the cheapest retained purchase
with the prescribed witness counts. This preserves the optimum
with the small witness proved above, and every accepted purchase
is feasible.

For \(r,k\ge1\), the number of failure factors is at most
\[
 H=\sum_{i=0}^{\min(k,m)}\binom mi\le (em/k)^k\le(2er)^k.
\]
The middle bound follows by evaluating \((1+x)^m\) at \(x=k/m\).
Each factor has total degree at most \(r(k+1)\), and there are
at most \(r(r+k)+m=O(r(r+k))\) shared variables. The dimension
of the coefficient space is at most
\((1+r(k+1)H)^{O(r(r+k))}\).
The logarithm of this quantity is \(O(kr(r+k)\log(r+1))\).
Polynomial multiplication and elimination take polynomial time in
this dimension, the factor count and the input size. Determinant
expansion uses at most \(r!\) terms, which the bound absorbs.
There are polynomially many target-count pairs and prefix positions,
and the new group state adds at most \(2(k+2)\) combinations per
witness state. Only the current group contributes counters.

For rational bit complexity, first retain a row basis and clear
column denominators. The resulting integer entries have bit length
polynomial in the input length by exact elimination. With integer tags
at most \(m\), the coefficient \(\ell_1\)-norm of a determinant factor
is at most \(r!\) times the product of the largest entry
\(\ell_1\)-norms in its columns. Its logarithm is polynomial in the
input length and in \(r,k\). Norms multiply submultiplicatively, so
the product coefficients have bit length at most \(H\) times that
bound. Exact elimination on their coefficient vectors takes polynomial
bit time in the dimension and coefficient bit lengths, by the
minor-size bounds for fraction-free elimination. The resulting
exponential factor is unchanged. A retained candidate carries its
actual purchase and witness.

At rank zero, choose the cheapest feasible option independently
for each group. For \(k=0\), the complement of a purchase must be
independent in the dual matroid and have at most \(|G_i|-1\)
elements of each group. This is weighted matroid intersection
with a partition matroid, maximizing the weight of the complement
[2]. It gives the polynomial-time edge case, including for the
dual state-matching oracle.

## Actuator reduction and one-sided sampling

Make the bipartite state graph with edges from left state \(i\)
to right state \(j\) for allowed state arcs \(i\to j\). Let
\(T\) be its transversal matroid on the right states, of rank \(\rho\),
and \(M=T^*\). The structural controllability criterion for a
surviving dedicated placement \(X\) is that its accessible states
cover the graph and that the state-input matching covers all
right states [3]. The matching condition holds exactly when
\(E\setminus X\) is matchable by state arcs: unmatched right
states in \(X\) can be covered by their own dedicated inputs.
Equivalently, \(X\) spans \(T^*\), of rank \(d=n-\rho\).
Accessibility holds exactly when \(X\) meets each source SCC.
These SCCs give the disjoint robust groups.

Assign independent indeterminates \(x_{ij}\) to the allowed entries
of the state matrix and zeros elsewhere. For every \(\rho\)-column
true primal basis \(D\), choose a matching that saturates \(D\)
and use its \(\rho\) rows. The resulting minor is a nonzero
polynomial of degree \(\rho\): its matching monomial is distinct
from all other permutation monomials. Conversely any nonzero
minor of a sampled matrix supplies a nonzero term and hence a
matching of its columns. The sampled rank is at most \(\rho\).

Pick each allowed entry uniformly from \(\{1,\ldots,2^b\}\) over
the rationals and discard a sample of rank less than \(\rho\).
At equal rank, every sampled primal basis is a true basis.
Represent its dual with rows forming a basis of the kernel of the
sampled matrix. For a \(d\)-set \(X\), projection of the kernel
onto the coordinates \(X\) is bijective exactly when the
complementary \(\rho\) columns are independent. This proves the
dual formula and shows that every sampled dual basis is a true
dual basis. Spanning in the sampled dual therefore implies
spanning in the true dual, also after every allowed failure.
The quota test uses actual purchases and statuses. A bad sample
can suppress an optimum, return a costlier feasible placement, or
return without a placement.

If \(\rho=0\), the dual is free and the only spanning purchase
contains all states. If \(d=0\), the dual has rank zero and groups
can be solved directly. Both cases skip sampling.

## Sampling for representation and for optimization

By Schwartz–Zippel [4], each prescribed matching minor has
probability at most \(\rho/2^b\) of vanishing. A union bound
over all true bases gives an error probability at most
\(\rho\binom{n}{\rho}/2^b\le \rho\,2^n/2^b\) for exact
representation. In particular,
\[
 b=n+\left\lceil\log_2(n/\delta)\right\rceil
\]
suffices for \(n\ge1\). The exponents in \(2^n\) and \(2^b\)
denote, respectively, the coarse bound on the number of bases and
the size of the sampling interval. On the event that all prescribed
minors survive, the rank check passes as well.

For optimization, fix a true minimum-cost purchase \(P^\star\)
and one of its small spanning witnesses \(C^\star\), independently
of the sample. For each allowed failure \(F\subseteq C^\star\),
choose a true dual basis in \(C^\star\setminus F\). Its
complement is a true primal basis with one prescribed matching
minor. If all these minors survive, the rank check passes and
\(C^\star\) still spans after every failure in the sampled
dual. Thus \(P^\star\) is feasible in the sample. Sampled
feasibility implies true feasibility, which forces equality of
optimal costs on this event.

The required bases can be reduced further. Within \(C^\star\),
take an inclusion-minimal family of bases such that each allowed
failure is disjoint from some family member. For each basis \(B_i\)
in the family there is a failure \(F_i\) avoided by \(B_i\) alone.
Thus \(B_i\cap F_i=\varnothing\), and \(B_i\cap F_j\ne\varnothing\)
for distinct \(i,j\). Write \(s=\min\{n,(k+1)d,|V|\}\) and
\(f=\min(k,s)\). The set-pairs inequality bounds the family size
by \(\binom{d+f}{f}\). Indeed, in a random ordering of \(C^\star\)
the events that every element of \(B_i\) precedes every element of
\(F_i\) are pairwise disjoint by the two cross-intersections.
Each event has probability at least \(1/\binom{d+f}{f}\).
The count is sharp for \(U_{d,d+k}\) with all elements vulnerable:
each \(k\)-failure leaves exactly one basis.

Hence one can use
\[
 A=\min\left\{\binom{n}{\rho},
              \sum_{i=0}^{f}\binom{s}{i},\
              \binom{d+f}{f}\right\}.
\]
The number of distinct required minors is at most \(A\). A single
trial, including rank rejection as an unsuccessful trial, has
optimality error at most \(\rho A/2^b\) on a feasible graph. For
\(\rho,d>0\), it suffices to use
\[
 b=\left\lceil\log_2\left(\frac{\rho A}{\delta}\right)\right\rceil
   =O\bigl(\log(n+1)+\log\binom{d+k}{k}+\log(1/\delta)\bigr).
\]
In particular \(b=O(\log(n+1)+d+k+\log(1/\delta))\).
The API sets \(\delta=2^{-t}\) for an integer \(t\ge1\) and
uses \(t+\lceil\log_2(\rho A)\rceil\) bits per entry.
These are unconditional probabilities. If resampling until rank
acceptance, acceptance probability is at least \(1-\rho/2^b\)
by any one fixed true basis minor. For \(2^b>\rho\) the error
conditional on acceptance is at most
\((\rho A/2^b)/(1-\rho/2^b)\).

Matching, SCC computation, and the kernel take polynomial time.
Kernel coefficients are ratios of integer minors with bit length
\(O(\rho(b+\log(n+1)))\) by the determinant expansion bound.
Applying the represented theorem at rank \(d\) proves the bit
complexity in the actuator theorem. An exact verifier can use the
purchased set for the SCC quotas and test state matching for every
allowed failure of the small witness.

For a deterministic alternative, the rational transversal representation
of Gurjar et al. [5, Theorem 1] costs
\(2^{O(\log^2(\rho+1))}\operatorname{poly}(n)\) time and
bit length. Dualizing and using the same optimizer gives that extra
factor in the deterministic running time. Obtaining a deterministic
bound fixed-parameter tractable in \(d+k\) alone remains a research
target for this route, since \(\rho\) can grow at fixed \(d\).

## Reproduction and references

The [rational coefficient DP](scripts/linear_accessibility_quotas.py)
uses the same transitions at \(k=0\). The
[graph solver](scripts/randomized_actuator_placement.py) exposes
rank rejection and an unsuccessful optimization as separate trial
statuses. Its `verify_placement` routine tests each returned witness
against state matchings of the original graph. Run `python scripts/audit_linear_accessibility_quotas.py`
from the repository root for exact enumeration and adversarial samples.

1. M. Kundu, S. Nandi, and S. Saurabh, [*When Elements Fail: Weighted Fault-Tolerant Matroid Bases*](https://arxiv.org/html/2611.16333v1), arXiv:2611.16333v1, 2026.
2. E. L. Lawler, [*Matroid intersection algorithms*](https://doi.org/10.1007/BF01681329), Mathematical Programming 9, 31–56, 1975.
3. S. Pequito, S. Kar, and A. P. Aguiar, [*A Framework for Structural Input/Output and Control Configuration Selection in Large-Scale Systems*](https://doi.org/10.1109/TAC.2015.2437525), IEEE Transactions on Automatic Control 61, 303–318, 2016.
4. J. T. Schwartz, [*Fast probabilistic algorithms for verification of polynomial identities*](https://doi.org/10.1145/322217.322225), Journal of the ACM 27, 701–717, 1980.
5. R. Gurjar, D. Lokshtanov, P. Misra, F. Panolan, S. Saurabh, and M. Zehavi, [*Quasipolynomial-Time Deterministic Kernelization and (Gammoid) Representation*](https://doi.org/10.4230/LIPIcs.MFCS.2025.54), MFCS 2025.
