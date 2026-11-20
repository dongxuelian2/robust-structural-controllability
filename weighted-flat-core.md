# Weighted fault-tolerant spanning with a cardinality cap and disjoint bases

Let \(M=(E,\mathcal I)\) be a rank-\(r\) matroid given by an independence oracle.
Each element has a nonnegative binary-encoded cost and is either safe or
vulnerable. A set \(B\) is resilient if it spans \(M\) after deleting any at
most \(k\) of its vulnerable elements. Write \(p=r+k\) and \(b=(k+1)r\).

Kundu, Nandi and Saurabh [1, Theorem 3.1] give a
\(2^{O(p^4)}\operatorname{poly}(n,L)\)-time algorithm for this model, where \(L\)
is the input bit length. Their flat-saturation construction gives the
unrestricted bound below. The additional parameters retain a cardinality cap
and exploit disjoint bases of the selected set.

## Bounds

For \(r,k\ge 1\), let \(r\le m\le b\), \(1\le s\le k+1\), and \(sr\le m\).
Consider the minimum-cost resilient set among sets of size at most \(m\)
that contain \(s\) pairwise disjoint bases of \(M\). Define

\[
T=\binom{r+k-1}{k},\qquad j=\min\{r,\lfloor(k+1)/s\rfloor\},\qquad
t=\begin{cases}
T,&s=1,\\
1+\min\bigl\{T,\ s\displaystyle\sum_{i=1}^j
        \binom{r}{i}\binom{m-r}{i-1}\bigr\},&s>1.
\end{cases}
\]

**Theorem.** A deterministic algorithm finds a minimum-cost set under these
restrictions or reports infeasibility. With
\(D=\log(m+1)+r\log(t+1)\), its total time is
\(2^{O(mD)}\operatorname{poly}(n,L)\) and its number of independence
queries is \(2^{O(rD)}\operatorname{poly}(n)\).

- For \(s=1,m=b\), the total time is
  \(2^{O(p^4)}\operatorname{poly}(n,L)\) and the oracle-query bound is
  \(2^{O(p^3)}\operatorname{poly}(n)\). The latter counts cached small-set
  queries throughout the search.
- For \(s=1,m=O(p^{3/2})\), the total time is
  \(2^{O(p^{7/2})}\operatorname{poly}(n,L)\).
- For \(m=b,s\ge\sqrt p\), we have
  \(\log(t+1)=O((1+(k+1)/s)\log(p+1))\). The total time is
  \(2^{O(p^{7/2}\log(p+1))}\operatorname{poly}(n,L)\).

The last two runs find an unrestricted optimum when one exists under their
respective restrictions. The caller supplies the cap or packing requirement.
Rank zero has the empty optimum. At \(k=0\), the allowed packing requirement
is one and matroid greedy finds a minimum-cost basis when \(m\ge r\).

## Proof

**Cardinality bound.** Suppose a resilient \(B\) contains \(s\) disjoint bases.
Use those bases as the first \(s\) layers, then successively take maximal
independent subsets of the remaining elements until there are \(k+1\)
layers. Their union \(B_0\) has size at most \(b\). A failure of at most
\(k\) elements misses some layer. Every omitted element lies in the closure
of that layer, so \(B\setminus F\subseteq\operatorname{cl}(B_0\setminus F)\).
Thus \(B_0\) is resilient, retains the designated bases, and costs at most
\(B\). This also proves the unrestricted cardinality bound with \(s=1\).
Loops can be removed from a feasible set while retaining resilience and all
designated bases.

**Obstructions to exchange.** Consider a feasible \(B\) of size at most
\(m\) with designated disjoint bases \(I_1,\ldots,I_s\). For \(x\in B\),
put \(Y=B\setminus\{x\}\). A failure \(F\) in \(Y\) for which
\(Y\setminus F\) falls short of rank \(r\) leaves rank \(r-1\), since adding
\(x\) restores rank \(r\). The hyperplane \(P=\operatorname{cl}(Y\setminus F)\)
avoids \(x\), is spanned by \(B\cap P\), and satisfies
\(|B\setminus P|\le k+1\).
There are at most \(T\) such hyperplanes by the
[set-pairs argument](weighted-fault-tolerant-fpt.md#sharp-deficient-hyperplane-count).

By disjointness, some \(I_\ell\) has at most \(j\) elements outside \(P\).
If exactly \(i\) elements are outside, \(1\le i\le j\), then \(I_\ell\cap P\)
has \(r-i\) elements. Extend it to a basis of \(P\) using \(i-1\) elements of
\((B\setminus I_\ell)\cap P\), which is possible because \(B\cap P\) spans
\(P\). Specifying \(\ell\), the omitted basis elements and the extension
determines \(P\). This gives the second hyperplane count in the definition
of \(t\). When \(s>1\) and \(x\) belongs to \(I_\ell\), also forbid
\(\operatorname{cl}(I_\ell\setminus\{x\})\). There is at most one such
additional hyperplane, and it avoids \(x\).

A same-status \(y\notin B\) outside all the forbidden hyperplanes yields
another feasible set \(Y\cup\{y\}\). When \(y\) survives, it repairs each
deficient survivor of \(Y\). When \(y\) fails, the failure with \(y\)
replaced by \(x\) in \(B\) has the same survivors. For \(s>1\), a
designated basis that contained \(x\) stays a basis after the swap, and
every other designated basis is unchanged. For \(s=1\), resilience itself
implies the existence of a basis.

**Weighted flat saturation.** Use the rule of [1, Lemma 5] separately on
the safe and vulnerable nonloops. Set \(L_0=0\) and
\(L_h=tL_{h-1}+m\). In increasing cost order, keep an element unless a
flat \(Q\) containing it already contains \(L_{\operatorname{rank}(Q)}\)
retained elements. This maintains the upper bound
\(|W\cap Q|\le L_{\operatorname{rank}(Q)}\) for every flat in each status
class. If \(x\in B\) was discarded in a rank-\(h\) flat \(Q\), its
discard witness consists of \(L_h\) retained elements no costlier than
\(x\). Every forbidden hyperplane avoids \(x\) and intersects \(Q\) in a
flat of rank at most \(h-1\), so it covers at most \(L_{h-1}\) of the
witness elements. At least \(m\) elements avoid the whole family. At most
\(m-1\) belong to \(B\setminus\{x\}\). Use an unused element for the
exchange. Iterating places an optimum in the union core, preserving its
size and packing.

**Construction and search.** Each saturated flat is spanned by its
retained elements. Otherwise their closure has smaller rank and the same
retained count, exceeding its strictly smaller cap. Thus the saturation
test enumerates independent subsets of the current retained set of size
at most \(r\), computes closure membership by independence queries, and
checks the count and membership of the incoming element. The core has
size \(N\le 2m\sum_{i=0}^{r-1}t^i\), with \(\log(N+1)=O(D)\).
The scan costs \(n(N+1)^{O(r)}\operatorname{poly}(r)\) operations and
oracle calls. Threshold integers have \(O(D)\) bits.

Enumerate core subsets of size at most \(m\). Check resilience over
their allowed failures. For \(s>1\), check packing by successively
removing a basis, memoizing by the available set and the number of
bases still required. There are at most \(s\,2^m\) such states per
candidate and at most \(2^m\) basis choices at each state. Each check
uses \(2^{O(m)}\) rank operations. A greedy rank scan needs only
independence queries on sets of size at most \(r\). Caching the answer
for each queried subset of the final core bounds the number of new
oracle calls during the entire search by
\(\sum_{i=1}^r\binom{N}{i}\le(N+1)^r\). Cache lookup,
candidate generation, and weight comparisons are included in the time
bound. The exchange proof justifies returning the cheapest candidate
or reporting infeasibility.

## Limits of the two restrictions

Both restrictions can fail for a connected matroid. For \(a\ge3\), take
the rank-\((2a-1)\) truncation of the direct sum of \(U_{a,3a-2}\) and
\(a\) rank-one classes of \(a\) elements each. Set \(k=2a-1\) and give
every element unit cost and vulnerable status. Before truncation the
rank is \(2a\). Losing rank in the uniform block takes at least
\(2a-1\) deletions. Losing rank in a rank-one class takes \(a\) deletions.
Thus \(2a-1\) deletions reduce the pre-truncation rank by at most one.

If an element of a rank-one class is omitted, deleting the rest of its
class and one entire other rank-one class costs at most \(2a-1\)
failures and loses two ranks. If an element of the uniform block is
omitted, deleting \(2a-1\) further uniform elements leaves at most
\(a-2\) there. Hence the whole ground set is the unique resilient set.
Its size is \(a^2+3a-2\). Every basis uses at least \(a-1\) uniform
elements, so there are at most three disjoint bases. To see connectedness,
take elements in different summands and extend them to a basis of the
direct sum. This basis is a circuit after truncation, so all summands
belong to the same component.

**Open problem.** Obtain a deterministic
\(2^{o((r+k)^4)}\operatorname{poly}(n,L)\) total-time bound for all
independence-oracle instances.

## Reference and implementation

[1] Madhumita Kundu, Soumi Nandi, and Saket Saurabh,
[*When Elements Fail: Weighted Fault-Tolerant Matroid Bases*](https://arxiv.org/html/2611.16333v1),
arXiv:2611.16333v1, 12 November 2026.

The reference implementation is [weighted_flat_core.py](weighted_flat_core.py).
Its `solve_independence` entry point takes a Boolean independence oracle.
Run `python audit_weighted_flat_core.py` to reproduce the finite checks.
