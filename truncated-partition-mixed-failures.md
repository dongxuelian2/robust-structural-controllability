# Exact polynomial-time mixed-failure spanning on truncated unit-partition matroids

Research note, November 16, 2026. **Scope:** a promised matroid class, not a stronger worst-case bound for arbitrary independence-oracle matroids.

## Problem, promise, and prior results

Let a rank-\(r\) matroid \(M=(E,\mathcal I)\) be accessed by an independence oracle. Elements have nonnegative binary-encoded integer costs and are partitioned into safe elements \(S\) and vulnerable elements \(V\). Given any \(k\ge0\), find a minimum-cost set \(B\subseteq E\) such that
\[
\operatorname{rank}_M(B\setminus F)=r\quad\text{for every }F\subseteq B\cap V,\ |F|\le k.
\]

**Promise.** The nonloop elements have parallel classes \(C_1,\ldots,C_p\), and a set \(I\) is independent precisely if it contains no loops, \(|I|\le r\), and \(|I\cap C_i|\le1\) for every \(i\). Thus \(M\), after deleting loops, is the rank-\(r\) truncation of a unit-capacity partition matroid. This includes every rank-one or rank-two matroid and uniform matroids of all ranks. The algorithm assumes this promise; it does not recognize the property on arbitrary oracle inputs.

Bentert--Fomin--Golovach--Morelle (ESA 2025, Proposition 5) proved an \(O(n^2)\) algorithm on this class in the **all-vulnerable** weighted setting. Kundu--Nandi--Saurabh (arXiv:2611.16333v1, 2026, Theorem 1.1) proved a \(2^{O((r+k)^4)}n^{O(1)}\) general **mixed-status** independence-oracle bound. Neither earlier result is claimed here.

**Theorem (polynomial-time mixed-status extension).** For matroids satisfying the promise, weighted fault-tolerant spanning with arbitrary safe/vulnerable element types and arbitrary failure budget \(k\) is solvable exactly in \(O(n^2)\) independence-oracle queries plus \(O(n^4\log(n+1))\) integer arithmetic operations, with polynomial bit complexity. A simpler reference dynamic program takes \(O(n^5)\) operations. The algorithm reports infeasibility or returns an optimum and a witness. The result is specific to this class; it does not improve the quartic general-oracle exponent.

## Normalization and failure threshold

Loops can be removed since they do not contribute to any surviving rank and weights are nonnegative. Parallel classes are discovered by singleton and pair independence tests, and \(r\) by a greedy scan of one representative from each class. Rank zero is solved by the empty set.

Put \(k'=\min(k,|V|)\), \(K=k'+1\le n+1\). In every class, it suffices to choose either (a) the cheapest safe element only, (b) at most \(K\) cheapest vulnerable elements, or (c) nothing. This follows because a selected safe element makes every other selected element of its class redundant, \(K\) selected vulnerable copies cannot all fail, and parallel elements are interchangeable at no increase in cost.

For each class \(C_i\) define \(d_i=K\) if the safe choice is used, and otherwise let \(d_i\) be the selected vulnerable multiplicity, \(0\le d_i\le K\). A surviving set has rank \(r\) exactly when at least \(r\) distinct classes remain occupied. An adversary must empty at least \(m=p-r+1\) classes to make the remaining rank smaller than \(r\). Emptying a vulnerable-only class costs \(d_i\) failures; a selected safe class cannot be emptied with \(K-1\) failures and can be assigned cost \(K\). Empty classes cost zero. With \(d_{(1)}\le\cdots\le d_{(p)}\), therefore,
\[
\boxed{\ B\text{ is resilient}\quad\Longleftrightarrow\quad
        \sum_{i=1}^{p-r+1}d_{(i)}\ge K.\ }\tag{1}
\]
This includes initially nonspanning sets, \(k=0\), and unavailable element copies. The condition is exact because the adversary's cheapest strategy is to empty the \(m\) classes with the smallest emptying costs.

## Exact order-statistic duality

For any integers \(0\le d_i\le K\), with \(m=p-r+1\),
\[
\boxed{\ \sum_{i=1}^{m}d_{(i)}
 =\max_{0\le t\le K,\ t\in\mathbb Z}
 \left[\sum_{i=1}^p\min(d_i,t)-(r-1)t\right].\ }\tag{2}
\]
To prove the upper bound, fix the indices of the smallest \(m\) entries. Their truncated contributions are at most their actual values, and the other \(p-m=r-1\) contributions are each at most \(t\). For equality, select any integer \(d_{(m)}\le t\le d_{(m+1)}\) (or \(t\ge d_{(p)}\) if \(m=p\)). In that interval exactly the desired \(m\) entries contribute their values and the others contribute \(t\). All data are integral and bounded by \(K\).

As \(K\ge1\), (1)--(2) imply:
\[
B\text{ is resilient}\ \Longleftrightarrow\
\exists t\in\{1,\ldots,K\}:\ 
\sum_{i=1}^p\min(d_i,t)\ge Q_t:=(r-1)t+K.\tag{3}
\]

## Exact DP and complexity

Fix \(t\). For class \(C_i\), the admissible (coverage, cost, chosen elements) alternatives are:

- \((0,0,\emptyset)\);
- \((t,w(s_i),\{s_i\})\), where \(s_i\) is its cheapest safe element, if any;
- \((j,\sum_{\ell=1}^j w(v_{i,\ell}),\{v_{i,1},\ldots,v_{i,j}\})\) for \(1\le j\le\min(t,|C_i\cap V|)\), where vulnerable elements are sorted by cost.

Selecting more than \(t\) vulnerable copies cannot improve the capped coverage for fixed \(t\) and only adds nonnegative cost. The class alternatives are mutually exclusive. Define \(D_i[z]\) as the minimum cost for the first \(i\) classes with total coverage capped at \(Q_t\) equal to \(z\). Initialize \(D_0[0]=0\) and \(D_0[z>0]=+\infty\). For every alternative with coverage \(a\) and cost \(c\), update
\[
D_i[\min(Q_t,z+a)]
 \leftarrow\min\{D_i[\min(Q_t,z+a)],D_{i-1}[z]+c\}.
\]
The cheapest \(D_p[Q_t]\), minimized across all \(t=1,\ldots,K\), is the optimum; if every such value is infinite, the instance is infeasible. Back-pointers reconstruct a selected set. Every chosen set satisfies (3) and hence (1); every feasible optimum admits the normalization and satisfies (3) at some \(t\). Thus the algorithm is exact.

For each \(t\), the DP uses \(O(pQ_t(t+2))\) arithmetic operations. Summing gives
\[
O\left(\sum_{t=1}^{K}pQ_t(t+2)\right)
 =O\bigl(p(r+1)K^3\bigr)\subseteq O(n^5).
\]
Singleton and representative-pair tests take \(O(n^2)\) oracle calls and rank determination takes \(O(n)\) more. The numeric value of \(k\) can exceed \(n\) without increasing the bound because \(k'\le n\). Costs are sums of at most \(n\) input weights and have polynomial encoding length. Reconstructing one threshold run uses \(O(p\max_tQ_t)\) memory; value-only runs use \(O(\max_tQ_t)\).

## Monge acceleration to O(n^4 log n)

The preceding DP is intentionally simple. We can improve its time bound without changing its optimal solutions. Fix a threshold \(t\) and write \(\widehat D_i[z]\) for the minimum cost, using the first \(i\) classes, of attaining **at least** \(z\) units of threshold coverage. Initialize \(\widehat D_0[0]=0\) and \(\widehat D_0[z>0]=+\infty\). This is equivalent to the capped-coverage DP for purposes of the final minimum \(\widehat D_p[Q_t]\).

For class \(i\), let \(m_i=\min(t,|C_i\cap V|)\). Write its vulnerable prefix costs as \(P_i(j)=\sum_{\ell=1}^j w(v_{i,\ell})\) for \(0\le j\le m_i\), with \(P_i(0)=0\). Because the individual vulnerable weights are sorted in nondecreasing order, \(P_i\) is *discrete convex*: \(P_i(j-1)+P_i(j+1)\ge2P_i(j)\). The vulnerable-only transition is
\[
U_i[z]=\min_{0\le j\le\min(z,m_i)}
      \left(\widehat D_{i-1}[z-j]+P_i(j)\right).\tag{4}
\]
Choices with \(j>z\) are no cheaper than \(j=z\) because the prefix costs are nondecreasing. Let \(u=z-j\); the transition matrix is
\[
A(z,u)=\widehat D_{i-1}[u]+P_i(z-u),
\qquad \max(0,z-m_i)\le u\le z.
\]
For \(z_1<z_2\) and \(u_1<u_2\), whenever the four entries are defined, discrete convexity gives
\[
A(z_1,u_1)+A(z_2,u_2)\le
A(z_1,u_2)+A(z_2,u_1).\tag{5}
\]
Thus the matrix is Monge. Its permissible column intervals have nondecreasing endpoints as the row \(z\) increases. If two row minimizers inverted their column order, their cross entries would both be permissible; inequality (5) would contradict one of the minima, choosing the leftmost minimizer in case of equality. Hence the **leftmost minimizing column is nondecreasing** across rows.

The standard divide-and-conquer monotone-row-minima search evaluates each of the \(Q_t+1\) rows using \(O(Q_t\log(Q_t+1))\) total matrix evaluations; unlike plain \(O(Q_t m_i)\) convolution, this bound does not depend multiplicatively on \(m_i\). The optional safe-element branch adds a single comparison per row:
\[
\widehat D_i[z]=\min\left\{U_i[z],\
 w(s_i)+\widehat D_{i-1}[\max(0,z-t)]\right\}.
\]
Omit the second branch if no safe element exists. Use a finite sentinel larger than the sum of all input weights to implement \(+\infty\), keeping the Monge ordering valid even for unreachable states. Store minimizing columns and class choices for witness reconstruction.

The accelerated algorithm uses
\[
O\left(\sum_{t=1}^K p Q_t\log(Q_t+1)\right)
=O\bigl(p(r+1)K^2\log(n+1)\bigr)
\subseteq O(n^4\log(n+1))
\]
integer arithmetic operations, still with \(O(n^2)\) oracle queries and polynomial bit complexity. This is a genuine reduction in the polynomial exponent **within the promised subclass**. It leaves the arbitrary-matroid oracle FPT exponent unchanged.

The file **truncated_partition_fast.py** implements this acceleration independently of the simpler **truncated_partition_fault.py**. Running **python audit_fast.py** checks 1,185 additional seeded fixtures against the baseline DP, including survivor ranks for every admissible failure of each returned witness.

## Example

For rank \(r=2\), failure budget \(k=2\), consider three classes. \(C_1\) has a safe element of cost 3 and a vulnerable element of cost 1; \(C_2\) has three vulnerable elements of cost 1 each; \(C_3\) has three vulnerable elements of cost 2 each. Selecting the safe element of \(C_1\) and all three vulnerable elements of \(C_2\) costs 6 and remains spanning after any two allowed failures. Without selecting the safe element, the cheapest feasible choice uses one vulnerable from \(C_1\), two from \(C_2\), and two from \(C_3\), costing 7. Hence the mixed-status optimum is 6.

## Reproducible verification

Run **python audit_truncated_partition.py** (Python 3.10+, standard library only). Its independent brute-force oracle enumerates every selected subset and every allowed failure for 462 seeded small fixtures, also checking the order-statistic identity; **303,728** selected subsets and all corresponding optimum values matched. Cases cover loops, zero weights, safe/vulnerable mixtures, all ranks from zero, and failure budgets exceeding the number of elements.

Finite computation supports but does not establish the theorem. Publication priority of this particular mixed-status extension requires additional literature review. The general-oracle \(2^{O((r+k)^4)}\) bound is **not** improved by this result. The solver uses only an independence oracle; its correctness requires the stated class promise, which is not verified for arbitrary oracle inputs.

## References

1. M. Bentert, F. V. Fomin, P. A. Golovach, L. Morelle, *Fault-Tolerant Matroid Bases*, ESA 2025, Proposition 5. https://doi.org/10.4230/LIPIcs.ESA.2025.83
2. M. Kundu, S. Nandi, S. Saurabh, *When Elements Fail: Weighted Fault-Tolerant Matroid Bases*, arXiv:2611.16333v1, November 12, 2026, Theorem 1.1. https://arxiv.org/abs/2611.16333
