# Clause constraints with rank-three matching blocks

Deciding whether the optimum equals twice the number of source SCCs is
NP-complete even when every component of the dual state matching matroid
has rank at most three. The graphs can be weakly connected, loopless and
bipartite, with underlying degrees two or three and no mandatory location.
The unpadded source SCCs have four or six vertices. The same reduction
admits girth at least N^(1-delta) for every fixed 0<delta<1, with the source
cycle lengths multiplied by a common odd factor.

The construction uses a seven-element cotransversal matroid that tests
one OR clause. The source quotas couple choices in different matching
blocks. For a formula with m clauses, maximum variable occurrence B and
minimum number u* of violated clauses, the cost satisfies

    gamma(B) = 1 + 2*floor((B-1)/2),
    7m + ceil(u*/gamma(B)) <= OPT <= 7m + u*.

These bounds prove NP-completeness at the accessibility lower bound and
transfer a bounded-occurrence gap. The resulting numeric hardness is
weaker than 250/249. The reduction also makes explicit where long variable
cycles lose a stronger source gap.

## The clause matroid

Use seven elements p,a1,b1,a2,b2,a3,b3. The dual T is presented by

    L1 = {p,a1,b1},    L2 = {p,a2,b2},    L3 = {p,a3,b3},
    L4 = {b1,b2,b3}.

The row rank is four. Each row has three elements, each pair has union
of size five, each triple has union of size at least six, and their
total union has size seven. After deleting a three-element set X, Hall's
condition for two or more rows still holds. A singleton row fails Hall
exactly if X is that row. The remaining three rows still match in that
case. The dual rank formula therefore gives rank three for every triple
of M=T* except the four L_i, which have rank two. Every pair is independent
and every four-set spans. Thus a four-set robustly spans M exactly when
it contains none of the L_i, and every set of size at least five is robust.

A robust four-set containing p must take one element from each {a_i,b_i}.
Take a_i to mean that literal i is true. The set containing all three b_i
is excluded by L4. The other seven patterns are exactly the satisfying
patterns of a three-literal OR clause. Adding one a_i repairs the all-false
pattern at cost one.

The three lines through p also certify nontransversality. They are cyclic
flats of rank two with pairwise unions of rank three and common
intersection of rank one. The Mason inequality would require 1<=0.
A rank-three transversal matroid on this ground set that retains all
original independent sets must make at least one of the lines independent.
Otherwise the lines would remain
cyclic flats with the same ranks and the same violation. If L_i is lost,
p+a_i+b_i+a_j for j different from i is a spurious robust four-set with
the anchor present. Thus fundamental assignment relaxations cannot
preserve this local clause test on all supports. See the cyclic-flat
criterion in [Bonin, Kung and de Mier, Theorem 3.2](https://arxiv.org/abs/1009.3435).

## Variable consistency and the cost bounds

Use one clause block per clause. Assume m is positive and even. For odd m
take two disjoint copies of the formula, including distinct copies of the
variables. This preserves B and the fraction of clauses that can be
satisfied. Put the m anchors into pairs, each with a quota of two.

For each literal occurrence o introduce a dummy d_o. Pair the 3m dummies
arbitrarily in the matching matroid, using a rank-one uniform block on each
pair. Robustness forces both members of every pair. Its dual presentation
has one row on the two dummies. The full matching matroid is the direct sum
of these blocks and the clause blocks.

Write t_o and f_o for the resources representing the two values of the
underlying variable. They are a_i,b_i in that order for a positive literal
and b_i,a_i for a negative one. Order the occurrences of each variable in
a cycle. For each occurrence o and its successor o' impose the quota

    |S intersect {t_o,f_(o'),d_o}| >= 2.

These groups and the anchor pairs partition the ground set. Since d_o is
selected, the constraint requires t_o or f_(o'). When both occurrences
are canonical, meaning exactly one of each pair is selected, the quota
forbids a transition from false to true. A cycle of canonical occurrences
therefore has a single Boolean value.

Every feasible set contains all m anchors and all 3m dummies. A clause
block costs at least four, so every placement costs at least 7m. At that
cost all clause blocks are canonical and each clause is satisfied.
Variable consistency gives a satisfying assignment to the formula.
Conversely, a truth assignment gives the canonical resources in each
clause and satisfies every quota. Repair each false clause with one
additional a_i. This preserves the quotas and gives cost 7m+u for an
assignment violating u clauses.

Consider a placement of size 7m+e. Call an occurrence double if it has
both resources and missing if it has neither. A clause block with local
excess h has d double and z missing occurrences, with d-z=h and d+z<=3.
If h=0 it is canonical. If h>0 these inequalities give d<=2h. There are
at most 2e double occurrences and at most e blocks with positive excess.

A variable with a missing occurrence must also have a double occurrence.
Indeed, without a double the occurrence after a missing one must be
canonical false. The canonical values cannot increase until the next
missing occurrence, which needs a preceding true resource. This would
violate its quota. Thus every variable without a double has consistent
canonical values. Round each variable by majority among its canonical
occurrences. A variable with a double has at most B-1 canonical
occurrences, so it disagrees with at most floor((B-1)/2) of them. At most
2e variables have doubles. A clause with zero excess was satisfied before
rounding, so it can become false only at a disagreement. Counting also
the at most e positive-excess clauses gives the bound gamma(B)*e.

The mandatory elements of this abstract instance are exactly the anchors
and dummies. Each literal element can be omitted from the full set while
retaining matching robustness and its group quota. The mandatory set has
rank 5m/2 out of the total matching rank 9m/2. Each clause supplies two
additional dimensions through its literal choices.

## Transfer to control and the quantitative scope

Apply `general-sparse-equivalence.md` to the cotransversal presentation
and the quota groups. The optimum is unchanged and the reverse map first
recovers a feasible abstract placement. The graph has 7m/2 source SCCs,
so 7m is also its accessibility lower bound. Its deficiency is 9m/2.
The clause blocks have rank three and the dummy blocks rank one. Row
splits act within blocks, and all other added matching elements are loops
or parallel copies. The target matching matroid still decomposes into
blocks of rank at most three. The two kinds of source groups yield cycles
of length four and six. There are no leaf completions in this family.
Before padding the graph has 113m/2-2 states and 143m/2-2 underlying edges.
Feasibility of a proposed placement is checkable by matching and SCC tests.

The bounded-occurrence source of Berman, Karpinski and Scott has B=4 and
m=1016s. Its stated gap has u*<=eta*s in the good case and u*>=(1-eta)*s
in the bad case for arbitrarily small fixed eta>0. Substitution gives
good cost at most (7112+eta)s and bad cost at least
(7112+(1-eta)/3)s. The limiting ratio is 21337/21336. This gives
constant-factor NP-hardness of approximation under the same graph and
matching-block restrictions using their source gap. See [Berman, Karpinski
and Scott, Theorems 1–2](https://theory.informatik.uni-bonn.de/ftp/reports/cs-reports/2003/85252-CS.pdf).

On formulas with three distinct variables in each clause, a random
assignment violates m/8 clauses in expectation. Conditional expectation
finds one with at most m/8 violations. The forward map therefore gives a
57/56 approximation on the constructed family when the formula is supplied.
This bound limits the OR construction as a route to a much larger gap.

## Why a long variable cycle loses the gap

For a literal ell, a bundle consists of the four clauses
ell OR (+/-s) OR (+/-t), with all four sign combinations and fresh
variables s,t. If ell is true all four are satisfied. If it is false,
every assignment violates exactly one clause of the bundle.

Take k bundles for NOT x, then a switch clause x OR y OR z with fresh y,z,
then k bundles for x, then one fresh satisfiable clause for even parity.
Use this occurrence order for x. There are m=8k+2 clauses with three
distinct variables per clause. Every global assignment violates at least
k clauses, and equality is attainable. Thus u*/m tends to 1/8.

Select the canonical x=false resources in the first phase and the
canonical x=true resources in the second phase. At the switch clause
select both x resources and one resource for each of y,z, plus the anchor.
This size-five set is robust. Both variable quotas adjacent to that
occurrence are satisfied. The transition back from true to false on the
wraparound is allowed by its quota. All remaining variables take a
consistent value, and all other clause blocks have size four. The cost
is exactly 7m+1, which is optimal because the formula is unsatisfiable.

Here B=8k+1. Consequently a single unit of excess can conceal a number
of violations linear in B. A stronger gap transfer from unbounded
occurrence formulas needs a different consistency construction.

`clause_quota.py` implements both truth maps. `verify_clause_quota.py`
checks all 128 clause ranks, all 31 bases, and 210 two-clause formulas.
It tests 176,610 candidate supports, rounds 85,658 feasible supports, and
lifts 800 truth assignments. Five sparse control fixtures check the
direct input-failure model, the ranks of all matching blocks and the
source quotas, including padded targets. Four sizes of the long-cycle
counterexample pass the source checks. Literature priority for the
combined control construction is unconfirmed.
