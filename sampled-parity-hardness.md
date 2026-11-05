# Hardness from sampling parity equations

**Theorem.** Unless $\mathrm{NP}=\mathrm{RP}$, for every fixed $\epsilon>0$
there is no polynomial-time $(92/91-\epsilon)$-approximation for unit-cost
1-RMSCP-d on general digraphs, or $(176/175-\epsilon)$-approximation on
strongly connected digraphs. Both results hold on loopless bipartite
graphs with underlying degrees two or three, no mandatory actuator
location and undirected girth at least $N^{1-\delta}$ for any fixed
$0<\delta<1$. The strongly connected targets are bidirected.

The reduction samples equations and works directly on the source
variables. The deterministic bounds from the amplified construction
are [843/838−ε and 1474/1469−ε](direct-lin2-hardness.md).

## Exact transfer of variable deletion

For a system $J$ of $s$ ternary parity equations on $n$ variables, write
$d(J)$ for the minimum number of variables to delete so that the
equations avoiding them are consistent. Use the shared backup and
private guards of the [direct parity reduction](direct-lin2-hardness.md),
with one packed job for each original variable and four guarded clause
jobs per equation. Each equation uses three distinct variables.

The selected labels contain at least one value at every variable.
Let $D$ be the variables with both labels selected. An equation touched
by $D$ has a selected literal in every one of its four clauses. For an
equation avoiding $D$, the single labels must satisfy the equation.
Conversely any deletion set and satisfying assignment on its complement
give exactly this pattern of selected labels. The shared backup repairs
the loss or diversion of one such label. The exact optima are therefore

$$
 \operatorname{OPT}_{\mathrm{gen}}=n+4s+2+d(J),\qquad
 \operatorname{OPT}_{\mathrm{sc}}=n+8s+2+d(J).
$$

Repeated equations can have separate guards and count separately in $s$.
Recover $D$ after the usual matching projection and dominated-resource
exchange. Those exchanges preserve selected original resources and
the paired quotas. The same realizations preserve the optima and graph
restrictions. The forced resources fall short of full matching rank
by $n-1$ in either construction.

## Sampling against all partial assignments

Use the regular Max-E3-Lin-2 source from Chlebík–Chlebíková, Definition 1
and Theorem 2. Each variable occurs equally often. For an arbitrarily
small fixed $\eta>0$, the optimum satisfied fraction is at least
$1-\eta$ in the good case and at most $1/2+\eta$ in the bad case.
Draw $s=3n$ equations independently and uniformly with replacement.
[Chlebík–Chlebíková, Theorem 2](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf).

In the good case fix a witness assignment before sampling. It violates
at most $3\eta n$ sampled equations in expectation. Deleting one variable
per violated equation suffices. For any fixed $\gamma>0$, choosing
$3\eta/\gamma<1/10$ gives $d(J)\le\gamma n$ with probability at least $9/10$.

In the bad case fix a set $D$ of $d=\lfloor n/7\rfloor$ variables and a
partial assignment to the remaining variables. Let $t$ be the source
fraction of equations touching $D$, and let $u$ be the source fraction
of equations that avoid $D$ and are violated. Regularity gives
$t\le3d/n$. Assign independent fair bits to $D$. Each touched equation
then fails with probability exactly $1/2$. Every full assignment
violates at least a $1/2-\eta$ fraction, so

$$
 u+\frac{t}{2}\ge\frac12-\eta,\qquad
 1-u\le\frac12+\eta+\frac{3d}{2n}\le\frac57+\eta.
$$

The chance this partial assignment survives the sample is at most
$(5/7+\eta)^{3n}$. Any smaller feasible deletion set extends to one
of size $d$. A union bound over the $\binom nd 2^{n-d}$ partial
assignments of this size gives probability at most

$$
 \exp\!\left(n\left[
 H(1/7)+\frac67\ln 2+3\ln(5/7+\eta)
 \right]\right),
$$

where $H$ is binary entropy with natural logarithms. The entropy term
uses that $H(x)+(1-x)\ln 2$ is increasing for $0\le x\le1/7$.
The exponent is negative for $\eta\le1/10000$, as certified by

$$
 \frac{7^7}{3^6}\left(\frac57+\frac1{10000}\right)^{21}<1.
$$

Thus $d(J)>n/7$ with probability tending to one. Substitution in the
exact optima gives limiting bad-to-good cost ratios $1+1/(7\cdot13)$
and $1+1/(7\cdot25)$ as $\gamma$ tends to zero. For a fixed approximation
factor strictly below either limit, choose $\gamma,\eta$ sufficiently
small and $n$ sufficiently large to distinguish the two cases by cost.
Finite smaller instances can be solved by exhaustive search.

The randomized reduction would put SAT in $\mathrm{BPP}$. Polynomially
many calls to the amplified SAT decision algorithm reconstruct a
satisfying assignment with high probability by self-reduction.
Verifying the reconstructed witness gives an $\mathrm{RP}$ algorithm
for SAT, which implies $\mathrm{NP}=\mathrm{RP}$. Uniform constraint
sampling for Max-CSP is prior work. The union bound here concerns
partial assignments and the exact cost of deleting variables.
[Trevisan, *Non-approximability results for optimization problems on bounded degree instances*](https://lucatrevisan.github.io/pubs/stoc01.pdf).

For a supplied source system, conditional expectation finds an assignment
violating at most $s/2$ equations. Deleting a variable per violation or
deleting every variable gives a selection with excess at most
$\min(n,s/2)$. Relative to the two baselines this gives factors at most
$10/9$ and $18/17$, respectively, over all source densities. These
limits apply to the guarded parity construction.

`verify_parity_guards.py` checks the exact transfer, the random-filling
identity and the rational exponent certificate. The $3/2-\epsilon$
target remains open.
