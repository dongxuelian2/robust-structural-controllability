# Hardness from stratified parity sampling

**Theorem.** Unless $\mathrm{NP}=\mathrm{RP}$, for every fixed $\epsilon>0$
there is no polynomial-time $(263/260-\epsilon)$-approximation for unit-cost
1-RMSCP-d on general digraphs, or $(503/500-\epsilon)$-approximation on
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

## Sampling three equations per variable

Use the regular Max-E3-Lin-2 source from Chlebík–Chlebíková, Definition 1
and Theorem 2. Each variable occurs exactly $\Delta$ times. For an arbitrarily
small fixed $\eta>0$, the optimum satisfied fraction is at least
$1-\eta$ in the good case and at most $1/2+\eta$ in the bad case.
Take $\Delta$ to be a sufficiently large power of two. Independently
for each variable, draw three uniform incident equations with replacement.
The sample $J$ has $s=3n$ equations counting multiplicity.
[Chlebík–Chlebíková, Theorem 2](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf).

In the good case fix a witness assignment before sampling. It violates
at most $3\eta n$ sampled equations in expectation, since every source
equation has the same expected number of draws. Deleting one variable
per violated equation suffices. For any fixed $\gamma>0$, choosing
$3\eta/\gamma<1/10$ gives $d(J)\le\gamma n$ with probability at least $9/10$.

## Bad-case bound with overlap among deletions

Put $p=3/20$ and fix a set $D$ of $d=\lfloor pn\rfloor$ deleted
variables and a partial assignment to the others. Write $x=d/n$.
Let $u$ be the source fraction of equations that avoid $D$ and are
violated, and let $t$ be the fraction touching $D$. If $m=n\Delta/3$
is the number of source equations, regularity gives

$$
 t=3x-c(D),\qquad c(D)=\frac1m\sum_e (|e\cap D|-1)_+.
$$

The sum counts repeated hits of the same equation. Filling $D$ with
independent fair bits violates each touched equation with probability
$1/2$. Soundness therefore gives $u+t/2\ge1/2-\eta$.
For variable $v$ let $u_v$ be the fraction of incident equations
that avoid $D$ and are violated. Then $\sum_v u_v=nu$ and
$u_v=0$ on $D$. By independence of the draws and the arithmetic–geometric
mean inequality, this partial assignment survives with probability at most

$$
 \prod_{v\notin D}(1-u_v)^3
 \le\left(1-\frac{u}{1-x}\right)^{3(n-d)}.
$$

Fix $\eta\le 1/10000$ and put $P=(1/2+1/10000+p/2)/(1-p)=5751/8500$.
The base in the last bound is at most $P-c(D)/(2(1-x))$.
Using $1-z\le e^{-z}$ and $1/(2P)>7/10>\ln2$ gives

$$
 \Pr(\text{survival})\le P^{3(n-d)}8^{-n c(D)}.
$$

To sum the overlap factor, now choose $D$ by independently deleting
each variable with probability $p$. Write $8^{-nc(D)}$ as the product
of the nonnegative factors $f_e=8^{-(n/m)(|e\cap D|-1)_+}$.
Each deletion bit occurs in $\Delta$ factors. Integrate the bits
one at a time, applying Hölder with exponent $\Delta$ to the factors
that depend on the current bit. Each such factor is replaced by its
conditional $L^\Delta$ norm. After all bits are integrated this yields
$\mathbb E\prod_e f_e\le\prod_e(\mathbb E f_e^\Delta)^{1/\Delta}$.
This is the product-space inequality of [Finner](https://doi.org/10.1214/aop/1176989534).
Since $\Delta n/m=3$ and each equation has three distinct variables,

$$
 \mathbb E\,8^{-nc(D)}\le G^{n/3},\qquad
 G=(1-p)^3+3p(1-p)^2+\frac{3p^2(1-p)}{512}+\frac{p^3}{512^2}.
$$

Any smaller feasible deletion set extends to one of size $d$.
Sum the survival bound over the $2^{n-d}$ assignments per set
of size $d$. Such a set has probability $p^d(1-p)^{n-d}$ in the
independent deletion measure. The union bound is at most

$$
 p^{-d}\left(\frac{2P^3}{1-p}\right)^{n-d} G^{n/3}
 \le\left[p^{-p}\left(\frac{2P^3}{1-p}\right)^{1-p}G^{1/3}\right]^n.
$$

The last inequality uses $d\le pn$ and $1/p>2P^3/(1-p)$.
The base is less than one by the exact rational certificate

$$
 \left(\frac{20}{3}\right)^9
 \left(\frac{40P^3}{17}\right)^{51}G^{20}<1 .
$$

Thus $d(J)>3n/20$ with probability tending to one. Substitution in the
exact optima gives limiting bad-to-good ratios $1+3/(20\cdot13)$
and $1+3/(20\cdot25)$ as $\gamma$ tends to zero. For a fixed approximation
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
and stratified survival bounds, and the rational exponent certificate.
The $3/2-\epsilon$ target remains open.
