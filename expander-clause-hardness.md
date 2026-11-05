# Clause gap transfer with bicircular variable copies

Let a formula consist of an even positive number $m$ of signed
three-literal constraints, all OR or all $\operatorname{OXR}(x,y,z)=x\lor(y\oplus z)$.
Write $u^*$ for the minimum number violated by an assignment. An odd
number of constraints can be handled by taking two disjoint copies of the
formula. The construction places each literal occurrence in a separate
robust assignment block, with spare blocks where needed for consistency.

**Theorem.** For either $d=2$ or $d=4$ there is a polynomial reduction
to unit-cost 1-RMSCP-d with

$$
 L+\lceil u^*/C_d\rceil\le\operatorname{OPT}\le L+u^*,\qquad
 C_2=185,\quad C_4=19.
$$

Here $M$ is the number of copy blocks, $L=(2d+3)M+m$ for OR and
$L=(2d+3)M+3m$ for OXR. The number of blocks per variable with
$\ell$ occurrences is $\ell+O(\sqrt{\ell})$, uniformly in $\ell$.
A truth assignment violating $u$ constraints lifts at cost $L+u$.
A placement of cost $L+e$ yields an assignment violating at most
$C_de$ constraints. There is no bound on variable occurrence.

The state graph is weakly connected, loopless and bipartite, with every
underlying degree in $\{2,3\}$, maximum indegree and outdegree three, and
no mandatory actuator location. For any fixed $0<\delta<1$ it can have
undirected girth at least $N^{1-\delta}$ on $N$ states. Each component of
the dual state matching matroid has rank at most $2d+2$ for OR and
$2d+4$ for OXR. Here $d$ controls the copy consistency graph.

## Two minimum choices in a robust assignment block

For $p\in\{3,4,5,6\}$ use jobs $0,\ldots,2p-1$. Each edge of a
resource graph can serve either endpoint, and each loop can serve its
unique endpoint. The assignment matroid is the bicircular matroid of this
graph. Its resources are the cycle $C_{2p}$, the common edges $F$, and
two loops at jobs 0 and 2. In order, the two alternating matchings of the
cycle and the common edges are

$$
 A_i=(2i,2i+1),\qquad B_i=(2i+1,2i+2\bmod 2p),\quad 0\le i<p,
$$

$$
 F_i=(2i+1,2i+4)\quad(0\le i<p-2),\qquad
 F_{p-2}=(2p-3,2p-1).
$$

Put the two loops in one quota-two group. Every robust selection costs
at least $2p+1$. With exactly that many resources, the selected ordinary
edges and loops form one connected component. Indeed each component on
$k$ jobs must have at least $k+1$ resources to withstand a deletion.
Removing the two loops leaves a tree in which only jobs 0 and 2 can be
leaves. Any other leaf would lose its only resource after one deletion.
Thus the tree is a Hamilton path from 0 to 2. Conversely any such path
with the two loops is robust. After a path edge fails, each component
can assign its loop to its endpoint and its edges to the other jobs
away from that endpoint. After a loop fails the remaining loop does the
same for the whole path.

There are exactly two such paths in each graph. The paths in the first
and second columns use $A\cup F$ and $B\cup F$, respectively.

| Ports | True | False |
| --- | --- | --- |
| 3 | 0,1,4,5,3,2 | 0,5,3,4,1,2 |
| 4 | 0,1,4,5,7,6,3,2 | 0,7,5,6,3,4,1,2 |
| 5 | 0,1,4,5,8,9,7,6,3,2 | 0,9,7,8,5,6,3,4,1,2 |
| 6 | 0,1,4,5,8,9,11,10,7,6,3,2 | 0,11,9,10,7,8,5,6,3,4,1,2 |

To verify exhaustiveness, start the path at 0 with unused set
$U=\{0,\ldots,2p-1\}\setminus\{0,2\}$. At vertex $v$ branch on every neighbor $w\in U$
and continue at $w$ with $U\setminus\{w\}$. If $U$ is empty accept
exactly when $v$ is adjacent to 2. Thus a minimum block chooses the
same Boolean label in all $p$ pairs $(A_i,B_i)$. Both labels select
every $F_i$.

## Directed consistency with disjoint quotas

For a variable with $\ell$ occurrences, use a directed multigraph of
indegree and outdegree $d$ on $n\ge\ell$ copies, with one permutation of
the copies for each index $i=0,\ldots,d-1$. Designate distinct copies for
the $\ell$ occurrences. The digraph has at least $\phi_d |U|$ outgoing
and incoming arcs for every set $U$ of at most half its vertices, with
$\phi_4=2/9$ and $\phi_2=1/92$.

For $d=4$ take $n=s^2$, the smallest square at least $\ell$. On
$\mathbb Z_s^2$ use the four permutations

$$
 (x+2y,y),\quad (x+2y+1,y),\quad (x,y+2x),\quad (x,y+2x+1).
$$

The underlying undirected multigraph has degree eight, counting loops
twice. The Margulis–Gabber–Galil bound gives second adjacency eigenvalue
at most $5\sqrt2$. Applying the Laplacian bound to the centered
indicator of $U$ gives boundary size at least
$(8-5\sqrt2)|U|(1-|U|/n)\ge(9/20)|U|$. The two directed cut sizes are
equal by degree balance, so each is at least $2|U|/9$. The cut condition
is vacuous at $n=1$. [Gabber–Galil](https://doi.org/10.1016/0022-0000%2881%2990040-4),
[formalized bound for these permutations](https://isa-afp.org/browser_info/current/AFP/Expander_Graphs/Expander_Graphs_MGG.html).

For $d=2$ take the smallest $n=4s^2\ge\ell$ and replace each vertex of
the preceding square construction by a directed four-cycle. Port $i$
carries the old $i$-th incoming and outgoing arcs. The external arcs
form one permutation and the cycle arcs form the other. For $U$ of
size at most $n/2$, let $H$ and $P$ be its full and partial clouds.
The undirected boundary $b$ is at least $2|P|$ by the cycle arcs, and
at least $(9/20)|H|-8|P|$ by the old cut of $H$. It follows that
$|U|\le4|H|+3|P|\le(400/9+3/2)b<46b$. Directed cuts have size $b/2$.
For $\ell=1$ use a single vertex with two self-arcs. Both constructions
have $n=\ell+O(\sqrt{\ell})$.

A normal copy uses $p=d+1$ pairs. For every $i$-th arc $u\to v$
impose quota two on $\{A_i(u),B_i(v),F_i(u)\}$. If both copies have
minimum cost, the group forbids the transition from false to true.
Each $A_i$ appears at one outgoing arc, each $B_i$ at one incoming
arc and each $F_i$ at one outgoing arc. These groups are disjoint
even for self-arcs and parallel arcs. Pair $d$ is reserved for literals.

## Clause groups and cost

Give each constraint one dummy resource. Pair these $m$ dummies into
rank-one two-resource assignment blocks, each serving its own job.
Robustness forces all dummies. For an OR constraint, use a quota-two
group on its dummy and the three resources at pair $d$ corresponding
to making its three literals true. The constraint holds on minimum
copies precisely when the group is satisfied.

For OXR the first literal has a larger copy, $p=d+2$. Use the same
first group for $x\lor y\lor z$. The second group expresses
$x\lor\neg y\lor\neg z$ using pair $d+1$ at $x$, the opposite
resources of pair $d$ at $y,z$, and the larger copy's $F_d$ as
filler. This filler is selected in both minimum choices. The groups
are disjoint because each literal occurrence has its own copy and all
resources used by the second group were unused by the first or by
consistency groups. The conjunction is OXR, and a false OXR constraint
violates exactly one group under a common truth assignment.

The sum of the copy minima and dummy costs is $L$, with $m$ larger
copies costing two extra each for OXR. A truth assignment selects
one minimum choice at each copy and every dummy. All consistency
groups hold. For each violated constraint add one missing literal
resource from its deficient group. Monotonicity of robust assignment
and the disjoint groups give a feasible selection at cost $L+u$.

For a feasible selection of cost $L+e$, at most $e$ copies exceed
their minimum. For a variable, let $D$ be these copies, and let $U$
be the smaller of the two classes of minimum copies. If $U$ is false,
every arc leaving $U$ goes to $D$. If $U$ is true, every arc entering
$U$ comes from $D$. This follows from the forbidden transition. Thus
$\phi_d |U|\le d|D|$. Set the variable to the majority minimum label,
breaking ties arbitrarily. Across all variables at most
$(1+d/\phi_d)e$ copies exceed the minimum or disagree with the assignment.
If all contacts of a constraint are minimum copies with the assigned
label, its groups certify the constraint. Each such group needs a literal
resource since it has at most one filler. Every violated constraint
therefore uses one of the counted copies, and each copy is the contact
of at most one constraint. This proves the reverse loss with
$C_d=1+d/\phi_d$.

## Exact realization and numerical gap

To present a copy as cotransversal, add one private resource per job.
For each original resource $a$, use the transversal row consisting of
$a$ and the private resources for jobs it can serve. These rows have
full rank by their distinct $a$ entries. Their dual is the augmented
assignment matroid by the [assignment transfer](assignment-hardness-transfer.md).
The original copy is robust, and every private singleton is dominated
by an original resource. The same dominated-resource exchange removes
all added private resources without increasing cost and retains every
selected original resource. Hence it preserves all quota groups.
The dummy pairs each use one row on their two resources. The resulting
dual is the direct sum of the copy and dummy blocks and is coloop-free.

Apply the [exact quota realization](general-sparse-equivalence.md).
In reverse, recover a quota-feasible set, remove private resources by
the exchange, and round. Row splits take place within each block, and
all remaining extensions add only matching loops or parallel elements.
Thus the bound on component rank survives. Source SCCs have lengths
four, six and eight before padding, for anchors, arcs and clause groups.

Replicate every constraint a sufficiently large constant number of
times on the same variables. This preserves the violated fraction and
makes $M/m$ arbitrarily close to three from above. Håstad's Theorems 6.5
and 6.18 give perfect completeness against unsatisfied fractions
$1/8-\eta$ for OR and $1/4-\eta$ for OXR, for every fixed $\eta>0$.
Using $d=4$, the respective limiting baselines per constraint are 34
and 36. The reduction gives NP-hardness within $5169/5168-\epsilon$
for OR and $2737/2736-\epsilon$ for OXR, for every fixed
$\epsilon>0$. [Håstad, *Some optimal inapproximability results*](https://people.kth.se/~johanh/optimalinap.pdf).

The bound replaces the occurrence-dependent loss of the
[variable cycle](clause-quota-hardness.md) by a constant. The
[shared backup construction](shared-backup-hardness.md) gives the
stronger general numerical bound 913/908−ε.
`verify_amplified_clause_quota.py` checks the four Hamilton path lists,
the exact minimum anchored robust sets, dominated exchanges, small
directed cuts, and both truth maps. Two control fixtures also check
source quotas and direct input-failure matching, including odd padding.
