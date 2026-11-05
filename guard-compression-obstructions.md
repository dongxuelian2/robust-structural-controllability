# Limits of compressing shared-backup parity checks

The [guarded parity reduction](direct-lin2-hardness.md) has one job
per label pair and one per Boolean clause. Two ways of saving guards
within this construction fail. Optional private guards can make every
instance attain the matching lower bound. Combining tests of independent
equations in a single auxiliary block still requires four jobs per
ternary equation when the shared backup pins the labels.

## Optional guards can absorb all doubled labels

For each variable $v$, a choice job can use either label $v^0,v^1$
or the global backup $b$. A guard job can use $b$ or its partner $g$.
Each clause job can use its literals or its own private resource.
A private resource is mandatory through a paired quota in the general
construction, or through a private guard job and buddy in the strongly
connected construction. An optional private resource has neither
forcing constraint. Let $B$ count one label per variable, $b,g$,
every private resource, and every mandatory buddy.

There are $B-1$ jobs, so any robust selection costs at least $B$.
In any feasible selection let $D$ be the doubled variables and $O$
the omitted private resources. Failure of $b$ forces a label in
every pair. It leaves exactly $|D|$ labels beyond the choice jobs,
and each job whose private resource is omitted needs one of these.
Consequently $|O|\le |D|$, and the selection costs $B+|D|-|O|$.

**Lemma.** Start with a full assignment $\sigma$ and double a set
$D$ of variables so that the selected labels cover all clauses.
Suppose each $v\in D$ has a distinct optional clause $H_v$ containing
both $v^{1-\sigma(v)}$ and some other variable's canonical label
$w^{\sigma(w)}$. There is a feasible selection of cost $B$.

Select all resources counted in $B$ except the private resources of
the $H_v$, and include the extra labels of $D$. For the nominal
matching use the canonical labels for the choices and the extra label
of $v$ for $H_v$. Use private resources for the remaining clauses
and guards. After loss of $b$ or $g$, assign its survivor to the
global guard. After loss of a canonical label, use $b$ for its choice.
After loss of an extra label, fill its $H_v$ with the other canonical
label and send $b$ to the displaced choice.

After loss of a selected clause resource or diversion to its
private guard, use a satisfying label for the clause. If that
label is canonical, send $b$ to its choice.
If it is extra, first repair the displaced $H_v$ with its other
canonical label and send $b$ to that choice. A private guard uses
the surviving member of its pair after either member fails. In the
quota construction the support meets all paired quotas. This covers
every single failure and uses one backup at most.

For a parity equation, index its clauses by the tuples that violate
it. A clause contains the labels opposite to its indexing tuple.

- **All parity guards optional.** For any parity system of arity at
  least three, take any full assignment and an inclusion-minimal set
  $D$ hitting its violated equations. Each $v\in D$ has a distinct
  violated equation with $v$ as its only member of $D$. Flip two
  other bits in that equation's canonical tuple. The result still
  violates parity, and its clause is served by $v$'s extra label and
  by the two flipped variables' canonical labels. These clauses
  serve as the $H_v$ of the lemma.
- **Two optional guards per ternary equation in the amplified
  construction.** Start with any assignment consistent across the
  copies. In each violated equation choose an optional clause indexed
  by a tuple different from the canonical tuple. Two distinct
  violating tuples differ on exactly two coordinates. Double the
  contact where they agree, and use that optional clause for its
  $H_v$. Each equation uses its own contacts, so the choices are
  disjoint. The argument allows all equality guards to be mandatory.
- **Optional equality guards in the amplified construction.** Keep
  all parity guards mandatory and start with any assignment consistent
  across the copies. Double one contact per violated equation. If
  both clauses on each consistency edge are optional, each doubled
  contact can use one of its incident consistency edges for a buffer.
  Match the doubled contacts to distinct incident edges by Hall's
  theorem. Each contact has degree two and each edge meets at most
  two doubled contacts. On the assigned edge exactly one clause
  contains the extra label of the doubled contact and the canonical
  label of its neighbor.

In each case the optimum is exactly $B$ for every source instance.
The exact [assignment transfers](direct-lin2-hardness.md#realization-and-the-constants)
preserve this value. The first case applies to the
[sampled construction](sampled-parity-hardness.md) when all its
clause resources are optional.

## A job bound with selectable auxiliary resources

Consider one local block with $h$ matching jobs and selectable
auxiliary resources confined to the block. Each external variable
supplies one canonical label with its own choice job. A label can
also serve local jobs. Only the shared backup can replace a
canonical label at its choice job. Thus every canonical label is
pinned to its choice when the backup fails. Outside the block there
are only the choice jobs and the global guard. Any auxiliary quotas
depend only on the auxiliary support. The block accepts an
assignment if some support within a fixed auxiliary budget is
robust. The backup may also be eligible for local jobs.

Let $A(\sigma)$ be the set of local jobs adjacent to at least one
canonical label under assignment $\sigma$. Acceptance is monotone
in $A$. To see this, fix a successful auxiliary support $S$. Failure
of the backup requires $S$ alone to match all $h$ jobs. Failure of
a canonical label is repaired by the backup while $S$ covers the
block. For failure of an element $a$ of $S$, the remaining $S-a$
can lose at most one matching unit. If it needs a label to recover
the last unit, let $j$ be the job covered by that label in a
successful matching. Any assignment with $j$ in its $A$ can use
one of its labels at $j$ and send the backup to that label's choice.
The other repairs use only the auxiliary resources and the backup.
The same support thus succeeds for every assignment with a superset
of $A(\sigma)$.

Each row indicator $j\in A(\sigma)$ is an OR of signed literals.
For each Boolean cube edge from an accepted assignment to a rejected
one, some row is true at the accepted end and false at the rejected
end. Otherwise monotonicity would accept both.

**One parity equation.** For parity on $k$ variables, there are
$k2^{k-1}$ such edges. A nonconstant row with $t<k$ literals is
false at $2^{k-t-1}$ rejected assignments and separates $t$ edges
at each, for a total of at most $2^{k-2}$. A full-width row
separates at most $k$ edges. Hence $h\ge4$ for $k=3$ and
$h\ge2k$ for $k\ge4$. Ordinary parity clauses attain these counts
for $k=3,4$.

**Independent ternary equations.** For $q$ equations on disjoint
triples, draw a uniformly random satisfying assignment. Each
accepted assignment has $3q$ rejected one-bit neighbors. A row
separates an edge at the accepted end exactly when it has exactly
one true literal, in which case it separates one edge. A constant
row separates none. Partial triples supply independent fair
literal bits. Each full triple independently contributes either
zero or two true literals, or one or three. A triple of the latter
type contributes exactly one with probability $3/4$. If there
is any such triple, the probability of exactly one true literal
in the row is at most $3/4$. If there are none, the full triples
must contribute zero and the fair bits have exactly one success
with probability at most $1/2$. Summing over rows gives
$3q\le 3h/4$, so $h\ge4q$, even when rows mix different triples.
Separate parity clauses attain the bound.

The job bound allows choices among auxiliary supports. It applies
when the loss of the backup pins all external labels, and the
conjunction must hold for every canonical assignment on the disjoint
triples. The [fixed-anchor bound](parity-gadget-obstruction.md)
uses a different hypothesis and also enforces presence of the label
pairs inside the parity gadget.

`verify_guard_compression.py` checks the rank-bound selections,
their transfers to control placements, and the row counts for
one and two disjoint triples.
