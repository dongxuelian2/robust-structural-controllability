# Weighted resilient spanning with disjoint robust accessibility groups

Let \(M\) be a rank-\(r\) independence-oracle matroid with nonnegative
binary-encoded costs. Let \(V\) denote its vulnerable elements, and let
\(S_1,\ldots,S_g\) be disjoint subsets of the ground set. Seek a
minimum-cost selection whose survivors span \(M\) and meet every group
after any deletion of at most \(k\) selected vulnerable elements.
Elements outside \(V\) are safe.

Put \(K=k+1\). A selected group is resilient exactly when it contains a
safe element or at least \(K\) vulnerable elements. For group \(i\),
let \(a_i\) be the cheapest safe cost, and let \(p_i\) be the sum of
its \(K\) cheapest vulnerable costs. Assign \(+\infty\) to an
unavailable option. A group with both options unavailable makes the
instance infeasible. Define \(u\) to count the groups with
\(a_i<p_i<+\infty\).

**Theorem.** For \(r\ge1\), write \(b=(k+1)r\) and
\(H=\binom{r+k-1}{k}\). The problem is solvable in

\[
2^{O(b(\log(b+u+1)+r\log(H+1)))}\operatorname{poly}(n,L)
\]

time, where \(L\) is the input bit length. In particular, when \(u=0\)
the exponential factor depends only on \(r+k\), with any number of
groups and arbitrary assignments of safe and vulnerable elements.
Rank zero requires only the cheapest feasible completion in each group.
For \(k=0\), the group constraint is the ordinary cardinality quota one
throughout, so the [weighted cardinality-quota theorem](weighted-disjoint-quotas.md)
also applies independently of \(u\).

## Completion costs

For a forced set \(C\), let \(g(C)\) be the minimum cost of a superset
satisfying all group requirements. In each group with no selected safe
element, compare the cheapest safe element with the cost of buying the
cheapest missing vulnerable elements needed to reach \(K\). A group
already meeting its requirement adds zero cost. These independent
completions include loops if they are useful for group survival.

Every feasible full selection contains a resilient spanning witness
\(C\) of size at most \(b\), by the [independent-layer
argument](weighted-fault-tolerant-fpt.md#a-small-optimum-exists-also-with-weights-and-safe-elements).
The optimum is therefore the minimum \(g(C)\) over resilient spanning
witnesses of size at most \(b\). Witnesses can exclude loops, and any
completion of a resilient witness preserves its survivor ranks.

## Prices for groups with a common marginal bound

Call a group *protected* if it is counted by \(u\). For a group with
\(p_i\le a_i\) and \(p_i<+\infty\), let \(\theta_i\) be its \(K\)-th
cheapest vulnerable cost. Assign adjusted price
\(\max(0,w(x)-\theta_i)\) to vulnerable \(x\), and \(w(x)-p_i\) to
safe \(x\). For a group with \(p_i=+\infty\) and \(a_i<+\infty\),
assign \(w(x)\) to vulnerable \(x\) and \(w(x)-a_i\) to safe \(x\).
An ungrouped element has price \(w(x)\). Write the adjusted price as
\(v(x)\). In each group and status the ordering by \(v(x),w(x)\)
refines the original weight ordering.

For every unprotected or ungrouped \(x\notin Y\),
\[
g(Y\cup\{x\})-g(Y)\ge v(x),
\]
with equality if \(Y\) avoids the group of \(x\), or if \(x\) is ungrouped.
To prove this, first suppose \(p_i\le a_i\) and \(p_i<+\infty\).
When \(Y\) lacks a safe element in that group, completing its
vulnerable quota costs at most \(p_i\): the original \(K\) cheapest
vulnerable elements supply enough missing elements. This completion
is at most the cheapest safe cost. Therefore forcing a vulnerable
element has the ordinary cardinality-quota marginal, with original
cutoff \(\theta_i\), and forcing a safe element saves at most \(p_i\)
in completion costs. Equality for an untouched group follows by
the cutoff formula or the baseline cost \(p_i\), respectively. If
\(Y\) already contains a safe element, forcing another element
adds its purchase cost. If vulnerable completion is unavailable,
a vulnerable element never reduces the outstanding cost of a safe
element, and forcing the first safe element saves exactly \(a_i\).
These cases prove the bound.

Within any group, including protected groups, replacing a forced
element by an unused same-status element of no greater weight does
not increase \(g\). Take a cheapest completion of the original
forced set. If the new element is already in the completion, it
remains feasible as it stands. Otherwise swap the elements. The
group counts of each status are preserved.

## Core and exchange

Use the [cardinality-quota saturation proof](weighted-disjoint-quotas.md#flat-saturation-core-with-color-caps)
on the unprotected groups and the ungrouped elements, separately
for both failure statuses. Interpret each group as a color and put
all ungrouped elements in one additional color. Process elements by
\((v(x),w(x),x)\), using the same per-color cap \(L_h=HL_{h-1}+b\)
and total cap \(A_h=HA_{h-1}+bL_h+b\), with \(L_0=A_0=0\).
All flat tests use the full ambient matroid.

Separately compress each protected group and status with the per-flat
cap \(L_h\) alone, in original weight order. The union of all these
cores has size at most \(2(A_r+uL_r)\).

Fix a resilient witness \(C\) of size at most \(b\), and \(x\) outside
the union core. Write \(Y=C\setminus\{x\}\). At most \(H\) deficient
survivor hyperplanes avoid \(x\). A same-status replacement outside
all these hyperplanes and outside \(Y\) retains resilience. For a
protected group or a per-color discard, the saturated flat contains
at least \(b\) retained same-group, same-status elements outside
these hyperplanes. Pick one outside \(Y\). It has no greater weight
and the same-group completion inequality applies. The same proof
for an ungrouped color uses additivity of ungrouped costs.

For a total-cap discard, exclude the ungrouped color and the colors
already touched by \(Y\). These are at most \(b\) colors, each
contributing at most \(L_h\) elements in the rank-\(h\) discard
witness. The hyperplanes cover at most \(HA_{h-1}\) more elements.
At least \(b\) elements remain in untouched unprotected groups.
Choose one such \(y\). Its completion marginal is \(v(y)\le v(x)\),
which is at most the marginal of forcing \(x\) into \(Y\).
Every exchange reduces the number of witness elements outside the
core while preserving resilience and the cost bound.

Enumerate core subsets of size at most \(b\), check survivor ranks
for all allowed failures, and compute the cheapest group completions.
As in the cardinality case, each saturated flat has an independent
basis among the retained elements for the relevant cap, so a scan
enumerates subsets of size at most \(r\) to find saturation witnesses.
The logarithm of the union-core size plus one is
\(O(\log(b+u+1)+r\log(H+1))\). The core construction and search
therefore satisfy the stated bound.

For dedicated actuator placement, take the dual state-matching
matroid, of rank equal to the matching deficiency, and use the source
strongly connected components as the groups. The matching and
accessibility tests apply to the same surviving selections. This
gives an exact algorithm parameterized by the matching deficiency,
failure budget, and \(u\).

## A declining marginal in a protected group

For \(K=2\), let a group contain a safe element of cost \(5\) and
two vulnerable elements \(x,y\) of cost \(4\) each. The four
completion values are \(g(\varnothing)=5\),
\(g(\{x\})=g(\{y\})=8\), and \(g(\{x,y\})=8\).
The singleton marginal of \(x\) is \(3\) in an untouched group
and \(0\) after forcing \(y\). A price equal to the untouched
marginal fails to lower-bound all marginals in this group.
The theorem handles these groups by preserving their identities.
The remaining parameterized problem asks for a time bound independent
of \(u\).

Run `python audit_mixed_accessibility_quotas.py` for 160 exact
optimum comparisons and 34,352 price checks, including protected-group
and total-cap discards.
