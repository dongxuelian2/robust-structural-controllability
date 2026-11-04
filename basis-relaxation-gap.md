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

## Finite checks and the nontransversal prototype

`python verify_basis_relaxation.py` checks every fundamental exchange in
this family for k=1,...,12 by independent assignment matching, verifies
the relaxed failure witnesses and the true rank deficit, and checks the
full rank/robust-optimum tables for k=1. The latter gives OPT 8 versus 6.

It also examines all 513 bases of the connected rank-five nontransversal
prototype from `cotransversal-hardness-transfer.md`. Its true optimum is
seven. The assignment relaxation has optimum six for 396 bases, including
B={0,1,2,3,6}, and optimum seven for the other 117. For that displayed bad
basis, selected mask 215 is relaxed-robust but loses true spanning after
deletion of element 0. Thus merely knowing that the matroid is
nontransversal does not determine the objective loss of a particular
basis relaxation, and the finite count is not a universal statement about
choosing a best basis.
