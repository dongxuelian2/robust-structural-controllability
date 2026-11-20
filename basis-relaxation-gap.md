# A fixed-basis assignment relaxation can lose a factor approaching 3/2

Let M be a coloop-free rank-d matroid and B a basis. Define the assignment
matroid N_B with one job per b in B. The basis resource b serves only its
own job. For e outside B, give it the jobs

    {b in B : B-b+e is a basis of M}.

These are the supports of its fundamental circuit, with e removed.
M-independent sets are N_B-independent. Indeed, each subset X of an
M-independent set lies in the closure of the union of its job
neighborhoods in B. Thus |X| is at most the size of that union, and Hall
gives the required matching. Both matroids have rank d, so
rho(N_B)<=rho(M). The converse independence implication need not hold:
preserving all exchanges from one basis does not preserve robust spanning
away from that basis.

The failure can be large even when M itself is a fundamental transversal
matroid, is simple and cosimple, and has no mandatory actuator. It is not
limited to nontransversal examples or to parallel resources.

**Proposition.** For every k>=1 there is such a matroid M and basis B with

    rank(M)=4k+1,       rho(M)=6k+2,       rho(N_B)=4k+2.

Thus rho(M)/rho(N_B)=(3k+1)/(2k+1) tends to 3/2. No cost-nonincreasing
reverse map from this relaxation to M can exist for arbitrary bases.
This is a relaxation obstruction, not a computational hardness result:
the displayed family has explicitly known optima.

## Construction and the true optimum

Set t=3k+1. Use k hub jobs and t local jobs v_i. For each i add two
resources a_i,b_i, each serving all hub jobs and v_i. First consider only
these 2t original resources. Every local job has exactly two neighbors,
so both must be selected. They form a robust set: after a resource in
group i fails, use its mate for v_i; in every other group use one resource
for its local job. The t-1>=k spare resources cover the k hubs. The
original robust optimum is therefore 2t.

Now add one private resource for each of the t+k jobs. By the
dominated-resource normalization in `assignment-hardness-transfer.md`,
the optimum remains 2t=6k+2. Call the resulting assignment matroid M.
Its private resources form a fundamental transversal basis.

M is simple. Each private resource serves a different singleton job, and
each nonprivate resource serves at least two jobs, so any two resources
have a two-job matching. M is cosimple as well: after deleting any two
resources, every job with a surviving private resource uses it. If one
or two private resources are missing, distinct original resources cover
those jobs. For two missing local privates their original groups differ;
for two hubs there are enough original resources; and a mixed pair has
the same immediate assignment. If a private and an original resource
are missing, at least one original neighbor of that job remains. Hence
the remaining full ground set always spans. In particular there is no
two-element cocircuit and no resource mandatory in every robust set.

## The bad basis and its exact relaxation

Take both a_i,b_i for i<=k and only a_i for i>k. These 2k+(t-k)=t+k
resources form B: use one resource of each full group for its local job
and the other for a different hub, and use each remaining a_i locally.

Let A consist of the 2k basis resources in the full groups. In N_B these
are 2k distinct anchor jobs. The other t-k=2k+1 basis resources a_i are
local jobs of the relaxation. Direct basis exchanges give:

- For i>k, both b_i and the physical private resource of v_i can serve
  every anchor job and the single relaxation job a_i.
- Every private hub resource and every private resource of v_i with
  i<=k can serve every anchor job.
- Basis resources have their defining singleton jobs.

For example, an outside resource in a group i>k replaces any anchor
resource by supplying the missing hub unit: its group then has one
resource for v_i and one for a hub. A private local resource frees its
group's original resource for that hub. It cannot replace a different
nonanchor group because that group's local job would lose its only
remaining neighbor in the basis. The other stated exchanges follow by
the same hub reassignment.

Select a_i,b_i for every i>k, a total of 4k+2 resources. After any failure,
use the surviving resource in its local group for that group's job and
one resource in every other group for its job. Each of the other 2k groups
has one surplus resource, and these serve the 2k anchor jobs. This is
robust in N_B. Rank 4k+1 gives the reverse lower bound 4k+2, proving the
claimed optimum.

This same selected set has rank only 3k+1 in M: it can reach the 2k+1
physical local jobs in its own groups and the k physical hub jobs. It
cannot reach any of the first k local jobs. Hence the relaxation can
even lose nominal spanning by k dimensions, despite recording every
single exchange from B correctly.

The exact cotransversal realization theorem applies separately to M and
N_B, including leaf completion and girth padding. It preserves their
values; it does not supply a map between these two different matroids.

## The same limit with resources serving at most two jobs

The large neighborhoods in the preceding construction are unnecessary.
For every k>=3 there is a simple, cosimple fundamental transversal M with

    rank(M)=4k,         rho(M)=6k,         rho(N_B)=4k+1,

where each resource serves at most two jobs and each job has at most
seven resource neighbors in the presentation of M. The ratio 6k/(4k+1)
again tends to 3/2. The fundamental presentation of N_B has larger
neighborhoods.

Start with a cycle on k hub vertices, take three copies of every edge,
and subdivide each copy once. Call the resulting simple bipartite graph
H. It has k hubs of degree six and 3k intermediate vertices of degree
two. Jobs are its 4k vertices. Original resources are its 6k edges, each
serving its two endpoints; add one private resource per job.

Every original edge resource is forced before adding the privates,
because every intermediate job has precisely its two incident edges.
The original full set is robust: after any one edge fails, H remains
connected and contains a cycle, so its remaining edges match all vertices.
The private-resource exchange therefore gives rho(M)=6k. Simplicity and
cosimplicity follow as before. Alternatively, the path-cover formula in
`variable-basis-path-cover.md` gives the same value: any path partition
has at most 2k edges, since every edge touches a hub and each hub has
path degree at most two. A partition into k three-vertex paths, each
joining two intermediate vertices through their common hub, and k
singleton intermediate vertices attains pc(H)=2k.

For B take the 2k edges of one subdivided cycle and one incident edge for
each of the other 2k intermediate vertices. It is a connected spanning
unicycle, hence a basis of the vertex-assignment matroid. Its 2k cycle
edges are the anchor jobs of N_B. Each outside edge or private resource
at a noncycle intermediate vertex has precisely all anchors and its own
branch edge as fundamental neighbors. Each private resource at a cycle
vertex has all anchors as neighbors. To verify these assertions, add the
resource to the unicycle: its unique bicircular circuit contains the full
cycle, together with its branch edge when there is one. Removing a
different branch edge would isolate that branch's job.

Select both original edges at every noncycle intermediate vertex and
one physical hub-private resource. This gives 4k+1 resources. In N_B each
branch uses one resource for its local job. If the hub-private fails,
the 2k branch surpluses cover all anchors; if a branch resource fails,
its group uses its survivor locally and the hub-private supplies the
missing anchor unit. This proves relaxed robustness, and rank gives the
matching lower bound. The same set has true rank only 3k: its resources
reach only the hubs and the noncycle intermediate vertices.

## Finite checks and the nontransversal prototype

`python verify_basis_relaxation.py` checks every fundamental exchange in
this family for k=1,...,12 by independent assignment matching, verifies
the relaxed failure witnesses and the true rank deficit, and checks the
full rank/robust-optimum tables for k=1. The latter gives OPT 8 versus 6.
It also verifies all 15,480 fundamental exchanges of the arity-two family
for k=3,...,12, its resource/job degrees, the explicit path partition,
the true rank deficit and every failure of the relaxed witnesses. The
smallest sparse case additionally checks every pair for simplicity and
cosimplicity.

It also examines all 513 bases of the connected rank-five nontransversal
prototype from `cotransversal-hardness-transfer.md`. Its true optimum is
seven. The assignment relaxation has optimum six for 396 bases, including
B={0,1,2,3,6}, and optimum seven for the other 117. For that displayed bad
basis, selected mask 215 is relaxed-robust but loses true spanning after
deletion of element 0. Thus merely knowing that the matroid is
nontransversal does not determine the objective loss of a particular
basis relaxation.

More precisely, the 117 exact bases are precisely those containing
element 5. Every robust set containing 5 has at least eight elements,
whereas the optimum is seven. Consequently every basis contained in any
optimal robust set has relaxed optimum six. Choosing a basis inside an
optimal solution therefore does not justify equality of the two optima.
The example still has other bases giving equality; it does not settle
whether every cotransversal matroid admits some exact basis relaxation.

## Maximizing column degrees does not select a reliable basis

The dense family also refutes a concrete polynomial basis-selection rule,
even when its dual presentation is maximal. Here maximal means that adding
any missing row-column incidence changes the transversal matroid.

For the dense family, a dual presentation has two identical rows

    F_i = {a_i,b_i, private(v_i)} union {all hub privates}

for each i. To see this, start with the usual private-basis dual rows
{a_i} union privates(N(a_i)) and {b_i} union privates(N(b_i)). Enlarging
each to F_i changes no independent set: within these two rows, any matched
a_i or b_i can be reassigned to its original row, while private matches
can use either row. Each F_i is a flat of M: it spans all hub jobs and
v_i, and every element outside F_i can additionally reach a different
local job.

Parallel-extend every element of the bad basis B to a class of
L=2t+1 copies; leave the other elements single. Write M' for this
extension and C_e for the copy class of e. A dual transversal presentation
of M' consists of the expanded old rows F_i and L-1 additional identical
rows C_e for every expanded class. Indeed, a matching saturating all new
rows consumes L-1 members of each class locally. At most one remains for
the old rows. Complementing a full matching therefore gives exactly the
bases consisting of one copy of each member of an old M-basis.

Every expanded column has degree L+1=2t+2. Every other column has degree
at most 2t. The d expanded classes are independent, so every basis
maximizing the sum of column degrees consists of one copy of each member
of B. This holds for every tie-breaking rule.

This presentation is maximal. Its old expanded rows are flats of M',
and each new class row is a rank-one flat, because M is simple. In
general a row R that is a flat of the dual M' cannot be enlarged: for
e outside R, put F=E(M')\R. Duality gives
r_T(F-e)=r_T(F), and a maximum matching on F-e uses no row R. Adding
the incidence (R,e) increases the matching rank on F by one.

The parallel-extension lemma gives rho(M')=rho(M)=6k+2. For any of the
maximum-degree bases, its fundamental assignment relaxation is N_B with
additional duplicate resources. A copy of a basis element has a singleton
neighborhood, and a copy of an outside element has the same fundamental
neighbors as that element. The dominated-resource lemma, applied to this
feasible assignment instance, preserves its optimum 4k+2. Thus even this
optimal degree-based basis choice in a maximal presentation loses a ratio
tending to 3/2. Parallel extension preserves cosimplicity and hence the
absence of mandatory elements, although M' is no longer simple.

The verifier checks these expanded presentations for k=1,2,3, including
their dual ranks on fixed and random supports, every fundamental exchange,
and the relaxed failure witnesses. It checks the base dual-presentation
identity on all 8,192 subsets for k=1. These are finite certificates for
the general argument above, not a hardness reduction.

## Every basis can admit spurious minimum-size placements

Choosing a better basis cannot in general make the whole robust feasible
family exact. This stronger family-preservation assertion is false even
for matroids that are both transversal and cotransversal. It is distinct
from asking whether some basis preserves the minimum objective value.

For r>=3 split 2r elements into disjoint sets A,D, each of size r. Give
an assignment presentation with r-2 common jobs, one A-only job and one
D-only job. Every A-resource serves all common jobs and the A-only job;
every D-resource serves all common jobs and the D-only job. Call the
resulting matroid M_r. Every r-subset meeting both groups is a basis;
the only dependent r-subsets are A and D, each a circuit-hyperplane of
rank r-1. All smaller subsets are independent. Complementing a basis
preserves this description, so M_r is self-dual. It is simple, cosimple
and connected, and is both transversal and cotransversal.

Take any basis B and put s=|B intersect A|, with 1<=s<=r-1. Its
fundamental assignment relaxation preserves the dependence of A exactly
when s=r-1. Otherwise every A-resource outside B has all r basis jobs
as neighbors, making A independent in N_B. The analogous assertion for
D requires r-s=r-1. These two equalities cannot both hold for r>=3.

Choose a group H whose dependence is lost, and any basis element b in
the other group. Then H+b is a robust set of size r+1 in N_B. To see this
directly, resources of H that lie in B have distinct singleton jobs,
the remaining H-resources serve all jobs, and b supplies its singleton
job after any H-resource fails. After b fails, H itself still matches
all r jobs. In M_r, deleting b leaves the nonspanning hyperplane H, so
the same selected set is not robust.

The true optimum is nevertheless r+1: choose r+1 elements with at least
two from each group. Every deletion leaves a mixed r-subset, hence a
basis. Rank gives the matching lower bound, so every N_B also has
optimum r+1. Thus every basis admits a spurious placement of optimum
cardinality, while its optimum value is exact. For r=3 this already gives
a six-element example and all 18 bases fail family preservation.

The verifier checks the complete rank tables and every basis for
r=3,4,5,6. A relaxation loses one circuit-hyperplane if B takes r-1
elements from one group, and both otherwise. Consequently it has exactly
r or 2r false robust subsets, all of size r+1; the verifier checks these
counts as well. This prevents using a fixed-basis relaxation as an exact
feasibility replacement under additional constraints without a separate
reverse-map argument. It gives no approximation-hardness constant and
does not answer the best-basis optimum-value question.
