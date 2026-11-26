# Derandomization and complexity boundaries for weighted dedicated actuator placement

Let \(G=(V,E_G)\) be a nonempty directed state graph on \(n=|V|\) states,
where each state \(v\in V\) has a nonnegative binary-encoded cost \(w(v)\)
and a failure status, either safe or vulnerable (\(V_{\text{vuln}}\subseteq V\)).
Let \(L\) denote the total input bit length, \(k\ge 0\) the failure budget,
and \(k'=\min(k,|V_{\text{vuln}}|)\).
In the bipartite state graph \(B_G=(V_L,V_R;A_G)\) with \(V_L=V_R=V\) and
\((u_L,v_R)\in A_G\) iff \((u,v)\in E_G\), let \(\rho\) be the maximum matching
size and \(d=n-\rho\) the matching deficiency.
Let \(G_1,\ldots,G_g\) be the source strongly connected components (SCCs) of
\(G\), and let \(u\) count the *protected* source SCCs: those satisfying
\(a_i<p_i<+\infty\), where \(a_i\) is the cheapest safe cost in \(G_i\) and
\(p_i\) is the sum of the \(k'+1\) cheapest vulnerable costs in \(G_i\) (taking
\(+\infty\) when an option is unavailable).

The [randomized actuator algorithm](linear-source-quotas.md) solves this
problem in \(2^{O(kd(d+k)\log(d+1))}\operatorname{poly}(n,L,\log(1/\delta))\)
bit operations with one-sided error \(\delta\), whereas direct deterministic
representation of \(B_G\) via Gurjar et al. [5] incurs a
\(2^{O(\log^2(\rho+1))}\) factor depending on the primal matching rank \(\rho\).
Meanwhile, the [oracle accessibility algorithm](mixed-accessibility-quotas.md)
avoids \(\rho\) at the expense of a \(u^{O(kd)}\) factor in \(u\).
This note eliminates the dependence on \(\rho\) when \(u\) is bounded,
decouples \(u\) from \(d\) and \(k\), gives \(\rho\)-free and \(u\)-free
deterministic algorithms for \(d\le 2\) and for primal transversal matroids,
and establishes exact oracle and strict-gammoid parity barriers for the
remaining regime.

## Main theorems

**Theorem 1 (\(\rho\)-free cut-covering compression).** Let \(d,k\ge 1\) and
put \(H=\binom{d+k'-1}{k'}\). An exact minimum-cost resilient dedicated
actuator placement on \(G\), or infeasibility, can be computed deterministically
in
\[
  2^{O\bigl(\log^2(\min(\rho,u)+1)+d^2\log^2(H+1)+kd(d+k)\log(d+1)\bigr)}
  \operatorname{poly}(n,L)
\]
bit operations. For \(d=0\) or \(k=0\) there are deterministic polynomial-time
algorithms. When \(u=0\) or \(u\le 2^{O(\sqrt{kd(d+k)\log(d+1)})}\), the time
bound is \(2^{O(d^2\log^2(H+1)+kd(d+k)\log(d+1))}\operatorname{poly}(n,L)\),
independent of \(\rho\). For arbitrary \(u\), the \(u\)-dependence improves from
\(u^{O(kd)}\) to \(2^{O(\log^2(\min(\rho,u)+1))}\) and decouples from
\(d\) and \(k\).

**Theorem 2 (Exact deterministic low-deficiency and primal algorithms).**
Weighted resilient dedicated actuator placement with arbitrary safe/vulnerable
elements, arbitrary nonnegative binary-encoded weights, and arbitrary
source-SCC counts \(u\) is solvable deterministically, independent of \(\rho\)
and \(u\):
1. For \(d=0\), in \(O(|E_G|\sqrt n+n\log n)\) arithmetic operations and
   polynomial bit time.
2. For \(d=1\), in \(O(|E_G|\sqrt n+n\log n+n(k'+1)^2)\subseteq O(n^3)\)
   arithmetic operations and polynomial bit time.
3. For \(d=2\), in \(2^{O(k^2)}\operatorname{poly}(n,L)\) bit operations.
4. For an explicitly presented rank-\(r\) primal transversal matroid \(T\)
   with disjoint accessibility groups and \(r,k\ge 1\), in
   \(2^{O(kr(r+k)\log(r+1))}\operatorname{poly}(n,L)\) bit operations
   (and polynomial time for \(r=0\) or \(k=0\)).

**Theorem 3 (Rank-oracle lower bound).** Fix integers \(k\ge 1\) and
\(r\ge 2\), put \(t=\lfloor r/2\rfloor\ge 1\), and let \(u\ge t+1\).
Any deterministic rank-oracle (or independence-oracle) algorithm for weighted
resilient spanning with \(u\) disjoint accessibility groups makes at least
\(\binom{u}{t}\) queries on rank-\(r\) matroids of size \(n=u(k+2)+k+(r\bmod 2)\)
with \(u\) protected groups. Any randomized algorithm succeeding with
probability at least \(2/3\) makes at least \(\frac{1}{3}\binom{u}{t}\) queries.
For \(r\ge 3\) and \(t=\lfloor (r-1)/2\rfloor\ge 1\), the same lower bounds
hold on loopless rank-\(r\) matroids of size
\(n=u(k+2)+k+1+((r-1)\bmod 2)\).

Theorem 3 rules out any \(\operatorname{FPT}(r,k)\) query bound (and any
\(f(r,k)\,n^{o(r)}\) query bound, since \(u=\Theta(n/k)\) yields
\(\Omega_{r,k}(n^{\lfloor r/2\rfloor})\) queries) in the rank-oracle model,
establishing that removing \(u\) from the parameter exponent requires
exploiting the linear or graph structure of \(T^*\). At rank two (with or
without loops), \(\binom{u}{1}=u\) is linear and \(O(n^2)\) rank queries
recover all loops and parallel classes, yielding the exact
\(\mathbb{Q}\)-representation and \(2^{O(k^2)}\operatorname{poly}(n,L)\)
algorithm of Theorem 2(3); superlinear query lower bounds \(\binom{u}{t}\)
with \(t\ge 2\) begin at \(r\ge 4\) (and \(r\ge 5\) for loopless matroids).

**Theorem 4 (Strict Gammoid \(q\)-Parity reduction).** Let \(q\ge 3\) be fixed.
Consider **Strict Gammoid \(q\)-Parity**: given a directed graph \(D=(V,E_D)\),
a sink set \(T\subseteq V\) of size \(d=qt\) (\(t\ge 1\)), and pairwise
disjoint \(q\)-element sets \(G_1,\ldots,G_g\subseteq V\), decide whether some
\(t\) groups have union forming a basis of the strict gammoid \(M(D,T)\).
This problem reduces in deterministic polynomial time to weighted dedicated
actuator placement on a state graph \(G^*\) with matching deficiency \(d\),
failure budget \(k=q-1\), \(u=g\) protected source SCCs, and polynomial integer
weights. Consequently, any deterministic \(\operatorname{FPT}(d,k)\) algorithm for
weighted dedicated actuator placement on general state digraphs yields a
deterministic \(\operatorname{FPT}(d)\) algorithm for Strict Gammoid
\(q\)-Parity (in particular, at \(k=2\), for Strict Gammoid \(3\)-Parity).

**Proposition 5 (Multilinear tensor obstruction in scale-\(t\) representations).**
Let \(\mathbb{F}\) be a field, \(V\) an \(\mathbb{F}\)-vector space of
dimension \(D=dt\), and \(1\le s\le (d-1)t\). Then
\[
  \operatorname{span}\{W^{\otimes s}: W\subseteq V,\ \dim W\le (d-1)t\}=V^{\otimes s}
  \quad\text{and}\quad
  \operatorname{span}\Bigl\{\bigwedge\nolimits^s W: W\subseteq V,\ \dim W\le (d-1)t\Bigr\}
  =\bigwedge\nolimits^s V.
\]

**Open problem.** Establish whether Strict Gammoid \(3\)-Parity admits a
deterministic \(\operatorname{FPT}(d)\) algorithm, and whether
weighted dedicated actuator placement admits a deterministic
\(\operatorname{FPT}(d,k)\) algorithm for unbounded \(u\) at \(d\ge 3\).

---

## Proof of Theorem 1: Ingleton–Piff equivalence and cut-covering compression

### Strict-gammoid presentation of the dual matching matroid

Let \(M=T(B_G)^*\) be the dual of the state transversal matroid on \(V\). By
[Pequito, Kar, and Aguiar](https://doi.org/10.1109/TAC.2015.2437525) [4, Theorem 3],
a surviving placement \(X\subseteq V\) is structurally controllable iff \(X\)
meets every source SCC \(G_1,\ldots,G_g\) and \(r_M(X)=d\). Fix a maximum
matching \(M_0\) of \(B_G\) of size \(\rho\). Let \(R_0\subseteq V\) and
\(L_0\subseteq V\) be the matched right and left sets, write
\(D_0=V\setminus R_0\) (\(|D_0|=d\)), and let \(\mu:R_0\to L_0\) be the
bijection satisfying \((\mu(r)_L,r_R)\in M_0\). Define the strict-gammoid
digraph \(D_H=(V,A_H)\) with out-arcs
\[
  (r,r')\in A_H \quad\Longleftrightarrow\quad r\in R_0,\ r'\ne r,\ (\mu(r)_L,r'_R)\in A_G.
\]
Vertices in \(D_0\) are sinks of \(D_H\). Recall that \(X\subseteq V\) is
independent in the strict gammoid \(M(D_H,D_0)\) iff \(X\) is linkable to
\(D_0\) by \(|X|\) vertex-disjoint directed paths in \(D_H\) (paths of length
zero allowed).

**Lemma 1 (Ingleton–Piff [1]).** \(M=M(D_H,D_0)\).

*Proof.* It suffices to verify equality of basis families, i.e., for every
\(B\subseteq V\) of size \(d\), \(V\setminus B\) is matched in \(B_G\) iff
\(B\) is linkable to \(D_0\) in \(D_H\).

Suppose first that \(\mathcal{P}=\{P_1,\ldots,P_d\}\) is a collection of
vertex-disjoint paths in \(D_H\) from \(B\) to \(D_0\). Because
\(|B|=|D_0|=d\), every vertex of \(D_0\) is a path endpoint and every vertex of
\(B\) is a path start. Let \(V(\mathcal{P})\) be the vertices on \(\mathcal{P}\).
Each \(r\in R_0\cap V(\mathcal{P})\) has a unique successor \(\sigma(r)\in V(\mathcal{P})\)
along its path, satisfying \((\mu(r)_L,\sigma(r)_R)\in A_G\). Define a matching
\(M'\subseteq A_G\) of size \(\rho\) by pairing \(\mu(r)_L\) with
\(\sigma(r)_R\) for \(r\in R_0\cap V(\mathcal{P})\) and pairing \(\mu(r)_L\)
with \(r_R\) for \(r\in R_0\setminus V(\mathcal{P})\). Every vertex of
\(V(\mathcal{P})\setminus B\) is the successor of some vertex in
\(R_0\cap V(\mathcal{P})\), and \(V(\mathcal{P})\supseteq D_0\), so the
right endpoints of \(M'\) are precisely
\((R_0\setminus V(\mathcal{P}))\cup(V(\mathcal{P})\setminus B)=V\setminus B\).
Thus \(B\) is a basis of \(M\).

Conversely, suppose \(B_G\) admits a matching \(M'\) of size \(\rho\) covering
\(V\setminus B\) on the right. Let \(L(M')\subseteq V_L\) be its left
endpoints. Consider the symmetric difference \(M_0\triangle M'\). Every
component of \(M_0\triangle M'\) is either an alternating cycle or an
alternating path. An odd-length alternating path would have both end-edges in
the same matching, making it augmenting for the other matching and contradicting
\(\nu(B_G)=|M_0|=|M'|=\rho\). Therefore every path component has even length and
has both endpoints on the same side of the bipartition. Toggling \(M'\) along
every path component whose endpoints lie in \(V_L\) preserves all right
endpoints \(V\setminus B\) (since internal \(V_R\)-vertices of those paths are
matched in both \(M_0\) and \(M'\)) and replaces \(L(M')\setminus L_0\) with
\(L_0\setminus L(M')\), producing a matching \(M''\) of size \(\rho\) with
right endpoints \(V\setminus B\) and left endpoints \(L_0\).

Let \(\phi:R_0\to V\setminus B\) be the bijection such that
\((\mu(r)_L,\phi(r)_R)\in M''\). Construct the directed subgraph
\(F=(V,\{(r,\phi(r)):r\in R_0,\ \phi(r)\ne r\})\subseteq D_H\). Every vertex in
\(F\) has outdegree at most \(1\) and indegree at most \(1\); vertices with
indegree \(0\) and outdegree \(1\) are precisely \(B\cap R_0\), vertices with
indegree \(1\) and outdegree \(0\) are precisely \(D_0\setminus B\), and
vertices in \(B\cap D_0\) are isolated. Hence \(F\) decomposes into directed
cycles and \(d\) vertex-disjoint directed paths from \(B\) to \(D_0\). Thus
\(B\) is a basis of \(M(D_H,D_0)\). \(_\blacksquare\)

### Core construction without \(u\) in the exponent

Assume \(d,k\ge 1\); the boundary cases \(d=0\) and \(k=0\) are resolved
below. Note first that \(k\) can be replaced by \(k'=\min(k,|V_{\text{vuln}}|)\)
without changing feasible sets or costs. If any source SCC has
\(a_i=p_i=+\infty\), return infeasibility. Otherwise set \(b=(k'+1)d\),
\(H=\binom{d+k'-1}{k'}\), \(L_0=A_0=0\), and
\[
  L_h = H L_{h-1}+b,\qquad A_h = H A_{h-1}+b L_h+b \qquad (1\le h\le d).
\]
Using the exact rank oracle \(r_M(Q)=|Q|-\rho+\nu_{B_G}(V\setminus Q)\) in
\(\operatorname{poly}(n)\) time per query, construct the flat-saturation union
core \(W\subseteq V\) of the [mixed-status accessibility theorem](mixed-accessibility-quotas.md):
1. On the union of all unprotected groups (\(p_i\le a_i\) or \(p_i=+\infty\)) and
   ungrouped elements, compress separately by status (safe and vulnerable)
   among nonloops in adjusted-price order \((v(e),w(e),e)\) with per-color cap
   \(L_h\) and total cap \(A_h\). Each status core has size at most \(A_d\), so
   enumerating independent subsets of size at most \(d\) within the current
   retained set takes \(O(n(A_d+1)^d)\) rank queries.
2. Compress each protected group \(G_i\) (\(1\le i\le u\)) **separately** by
   status among nonloops in weight order \((w(e),e)\) with per-flat cap \(L_h\)
   alone. Each group-status core has size at most \(L_d\), so processing all
   \(u\) protected groups takes \(O(u n (L_d+1)^d)\le O(n^2(A_d+1)^d)\) rank
   queries.

Thus \(W\) is constructed deterministically in
\(2^{O(d\log(A_d+1))}\operatorname{poly}(n)\) time, with no dependence on
\(\rho\) or \(u\) in the exponent, and satisfies \(|W|\le \min\{n,\,2(A_d+uL_d)\}\).
Check whether \(r_M(W)=d\); if \(r_M(W)<d\), the instance is infeasible by
[core optimum preservation](mixed-accessibility-quotas.md#core-and-exchange).

### Cut-covering compression to a smaller bipartite presentation

Form the terminal-disjoint augmentation \(\widetilde{D}_H\) on vertex set
\(V\cup S\cup T\), where \(S=\{w_{\text{in}}:w\in W\}\) and
\(T=\{t_{\text{out}}:t\in D_0\}\) are fresh disjoint copies of \(W\) and
\(D_0\), adding arcs \((w_{\text{in}},w)\) for \(w\in W\) and
\((t,t_{\text{out}})\) for \(t\in D_0\). Because vertices in \(S\) have
indegree \(0\) and unique out-neighbor \(w\), and vertices in \(T\) have
outdegree \(0\) and unique in-neighbor \(t\), a set \(X\subseteq W\) is
linkable to \(D_0\) in \(D_H\) iff \(X_{\text{in}}=\{w_{\text{in}}:w\in X\}\)
is linkable to \(T\) in \(\widetilde{D}_H\).

Apply the deterministic polynomial-time cut-covering algorithm of Koana and
Kumabe [8, Theorem 1.2] to \((\widetilde{D}_H,S,T)\), obtaining a vertex set
\(Z\subseteq V(\widetilde{D}_H)\) containing \(S\cup T\) of size
\[
  |Z|=O\bigl((|S|+|T|)^3\bigr)=O\bigl((|W|+d)^3\bigr)
\]
that contains a minimum \((A,B)\)-vertex cut in \(\widetilde{D}_H\) for every
\(A\subseteq S\) and \(B\subseteq T\). Define the torso digraph \(D_Z\) on
vertex set \(Z\), placing arc \((x,y)\) for distinct \(x,y\in Z\) iff
\(\widetilde{D}_H\) has a directed path from \(x\) to \(y\) whose internal
vertices all lie outside \(Z\). Vertices in \(T\) still have outdegree \(0\) in
\(D_Z\).

**Lemma 2 (Mimicking strict gammoid).** \(M(D_Z,T)|_S = M(\widetilde{D}_H,T)|_S \cong M|_W\).

*Proof.* Fix any \(A\subseteq S\). Any \((A,T)\)-cut \(C\subseteq Z\) in
\(D_Z\) hits every path from \(A\) to \(T\) in \(\widetilde{D}_H\): the
sequence of \(Z\)-vertices along such a path forms a walk from \(A\) to \(T\) in
\(D_Z\) because consecutive \(Z\)-vertices are joined by \(Z\)-internal-free
subpaths. Conversely, a minimum \((A,T)\)-cut \(C^*\subseteq Z\) of
\(\widetilde{D}_H\) hits every path from \(A\) to \(T\) in \(D_Z\), because
expanding each \(D_Z\)-arc along its \(Z\)-disjoint path in \(\widetilde{D}_H\)
cannot introduce vertices of \(C^*\subseteq Z\). By Menger's theorem, the
maximum number of vertex-disjoint paths from \(A\) to \(T\) is identical in
\(D_Z\) and \(\widetilde{D}_H\). \(_\blacksquare\)

Now construct the bipartite graph \(B_Z=(L_Z,R_Z;A_Z)\) with
\(L_Z=Z\setminus T\), \(R_Z=Z\), and
\[
  A_Z = \{(z_L,z_R):z\in Z\setminus T\}\cup\{(z_L,y_R):z\in Z\setminus T,\ y\in Z,\ (z,y)\in E(D_Z)\}.
\]
The matching \(\{(z_L,z_R):z\in Z\setminus T\}\) covers \(L_Z\), so \(B_Z\) has
maximum matching size \(\rho_Z=|Z|-d\), and applying Lemma 1 in reverse yields
\(T(B_Z)^*=M(D_Z,T)\). Restricting \(T(B_Z)^*\) to \(S\) therefore produces an
exact presentation of \(M|_W\).

### Representation, zero-padding, and complexity analysis

Compare \(\rho\) with \(\rho_Z=|Z|-d\): set \(B^*=B_G\) (rank \(\rho^*=\rho\),
columns \(m^*=n\)) when \(\rho\le\rho_Z\), and set \(B^*=B_Z\) (rank
\(\rho^*=\rho_Z\), columns \(m^*=|Z|\)) otherwise. In both cases
\(\rho^*=\min(\rho,\rho_Z)\) and \(m^*-\rho^*=d\).

Run the deterministic rational representation algorithm of Gurjar et al.
[5, Theorem 1] on \(B^*\). In
\((\rho^*)^{O(\log(\rho^*+1))}\operatorname{poly}(m^*)=2^{O(\log^2(\rho^*+1))}\operatorname{poly}(n)\)
time, it outputs a matrix \(A\in\mathbb{Q}^{\rho^*\times m^*}\) representing
\(T(B^*)\) whose entries have bit length \(B_0=2^{O(\log^2(\rho^*+1))}\).
Select \(\rho^*\) pivot columns of \(A\) and compute a kernel basis
\(K\in\mathbb{Q}^{d\times m^*}\) via Gaussian elimination. Clearing column
denominators in the \(\rho^*\times(\rho^*+1)\) submatrix for each free
coordinate yields integer matrix entries of bit length at most
\(\rho^*\log\rho^*+\rho^*(\rho^*+1)B_0=2^{O(\log^2(\rho^*+1))}\) by Hadamard's
inequality, in \(2^{O(\log^2(\rho^*+1))}\operatorname{poly}(n)\) bit
operations.

Assign column \(v_e\in\mathbb{Z}^d\) to each \(e\in W\) (from \(K\) directly
when \(B^*=B_G\), or from \(K[S]\) when \(B^*=B_Z\)), and assign the zero
column \(0\in\mathbb{Z}^d\) to every \(e\in V\setminus W\). Call the resulting
rank-\(d\) represented matroid \(M'\). Because \(M'|_W=M|_W\) and
\(r_M(W)=d\), every \(X\subseteq V\) satisfies
\[
  r_{M'}(X)=r_{M'}(X\cap W)=r_M(X\cap W)\le r_M(X).
\]
Consequently, any placement \(P\subseteq V\) feasible for \((M',G_1,\ldots,G_g)\)
also satisfies \(r_M(P\setminus F)=d\) and every source-SCC quota after every
allowed failure \(F\); hence \(\mathrm{OPT}(G)\le \mathrm{OPT}(M')\). Conversely,
by [core optimum preservation](mixed-accessibility-quotas.md#core-and-exchange),
\(W\) contains a resilient spanning witness \(C^\star\subseteq W\) whose
cheapest group completion \(P'\supseteq C^\star\) satisfies
\(w(P')=\mathrm{OPT}(G)\) and \(r_M(C^\star\setminus F)=d\) for all
allowed failures \(F\subseteq P'\cap V_{\text{vuln}}\). Because
\((P'\setminus F)\cap W\supseteq C^\star\setminus F\), we have
\(r_{M'}(P'\setminus F)=d\), so \(\mathrm{OPT}(M')=\mathrm{OPT}(G)\).

Execute the [represented-matroid DP](linear-source-quotas.md#dynamic-program-and-complexity)
on \((M',w,k',V_{\text{vuln}},G_1,\ldots,G_g)\). Its bit complexity is
\(2^{O(kd(d+k)\log(d+1))}\operatorname{poly}(n,L+L_{\text{rep}})\), where
\(L_{\text{rep}}=2^{O(\log^2(\rho^*+1))}\). From
\(\rho^*\le\min\bigl(\rho,\,O((A_d+uL_d+d)^3)\bigr)\le O\bigl((\min(\rho,u)+1)^3(A_d+L_d+d)^3\bigr)\)
and \(\log(A_d+1)=O(\log(k'+1)+d\log(H+1))\), we obtain
\[
  \log^2(\rho^*+1)
  =O\bigl(\log^2(\min(\rho,u)+1)+d^2\log^2(H+1)+\log^2(k'+1)\bigr).
\]
Both \(\log^2(k'+1)\) and the core-construction exponent \(O(d\log(A_d+1))\)
are absorbed into \(d^2\log^2(H+1)+kd(d+k)\log(d+1)\). This establishes the
bound in Theorem 1.

At \(k=0\) (or \(k'=0\)), a selection \(P\subseteq V\) is feasible iff
\(V\setminus P\) is independent in \(T(B_G)\) and omits at least one element
of each \(G_i\); maximizing \(w(V\setminus P)\) is weighted matroid
intersection of \(T(B_G)\) with a partition matroid, solvable in deterministic
polynomial time via Lawler [2]. The \(d=0\) case is proved below. \(_\blacksquare\)

---

## Proof of Theorem 2: Exact deterministic algorithms without \(\rho\) or \(u\)

### Deficiency zero

When \(d=n-\rho=0\), \(M=T(B_G)^*\) has rank \(0\) and every subset spans \(M\).
Only the pairwise disjoint source-SCC quotas remain. In each source SCC \(G_i\),
compare the cheapest safe state (if any) with the \(k'+1\) cheapest vulnerable
states (if at least \(k'+1\) exist). Returning the union of cheapest options,
or infeasibility if any \(G_i\) lacks both, takes \(O(|E_G|\sqrt n+n\log n)\)
arithmetic operations (maximum matching via Hopcroft–Karp, Tarjan SCCs, and
sorting).

### Deficiency one

When \(d=1\), \(D_0=\{t_0\}\) is a singleton. By Lemma 1, a state \(v\in V\)
is a nonloop of \(M=T(B_G)^*\) iff \(D_H\) contains a directed path from \(v\)
to \(t_0\). One reverse breadth-first search from \(t_0\) in \(D_H\) classifies
all nonloops in \(O(n+|E_G|)\) time after a single maximum matching of \(B_G\).
Every nonloop has rank \(1\) and all nonloops are mutually parallel. Therefore
a purchase \(P\subseteq V\) satisfies the spanning condition after every
allowed failure iff \(P\) contains at least one safe nonloop or at least
\(k'+1\) vulnerable nonloops.

Process the source SCCs \(G_1,\ldots,G_g\) and the ungrouped set
\(G_0=V\setminus\bigcup_i G_i\) in sequence. Maintain a dynamic-programming
table indexed by \(k'+3\) states: integer counts \(0,1,\ldots,k'+1\) recording
the number of selected vulnerable nonloops capped at \(k'+1\) (when no safe
nonloop has been selected), plus a single absorbing state \(\mathrm{SAFE}\)
recording that at least one safe nonloop has been selected.

Within a required group \(G_i\), sort its four disjoint subclasses
\(\mathrm{SN}_i\) (safe nonloops), \(\mathrm{SL}_i\) (safe loops),
\(\mathrm{VN}_i\) (vulnerable nonloops), and \(\mathrm{VL}_i\) (vulnerable
loops) by nondecreasing cost. By nonnegativity of weights and parallel
interchangeability:
- Choosing the cheapest safe nonloop \(\mathrm{SN}_i[0]\) satisfies \(G_i\) and
  moves every state to \(\mathrm{SAFE}\).
- Choosing the cheapest safe loop \(\mathrm{SL}_i[0]\) satisfies \(G_i\) and
  may be paired with the \(a\) cheapest elements of \(\mathrm{VN}_i\) for any
  \(0\le a\le \min(k'+1,|\mathrm{VN}_i|)\), shifting count \(v\) to
  \(\min(k'+1,v+a)\).
- Without a safe purchase in \(G_i\), group survival requires at least \(k'+1\)
  vulnerable purchases in \(G_i\); since the global vulnerable-nonloop counter
  is capped at \(k'+1\), choosing \(a\) cheapest elements of \(\mathrm{VN}_i\)
  and \((k'+1)-a\) cheapest elements of \(\mathrm{VL}_i\) suffices for each
  \(0\le a\le\min(k'+1,|\mathrm{VN}_i|)\) with
  \(k'+1-a\le|\mathrm{VL}_i|\).

For the ungrouped block \(G_0\), transitions consist of skipping, buying
\(\mathrm{SN}_0[0]\), or buying \(\mathrm{VN}_0[:a]\) for
\(1\le a\le\min(k'+1,|\mathrm{VN}_0|)\). At the end, take the minimum cost of
states \(k'+1\) and \(\mathrm{SAFE}\). Each group \(G_i\) has
\(O(\min(k'+1,|G_i|))\) transitions evaluated across \(k'+3\) states, taking
\(O((k'+1)\sum_i\min(k'+1,|G_i|)+(k'+1)^2)=O(n(k'+1)^2)\) operations. Together
with matching, SCCs, and sorting, the total arithmetic work is
\(O(|E_G|\sqrt n+n\log n+n(k'+1)^2)\subseteq O(n^3)\) with polynomial bit
lengths.

### Deficiency two

When \(d=2\), \(M=T(B_G)^*\) has rank \(2\). Identify all nonloops using the
dual rank formula \(r_M(\{e\})=1-\rho+\nu_{B_G}(V\setminus\{e\})\), i.e.,
\(\nu_{B_G}(V\setminus\{e\})=\rho\). Partition the nonloops into parallel
classes \(C_1,\ldots,C_m\) (\(m\le n\)) by testing
\(r_M(\{e,f\})=1\iff \nu_{B_G}(V\setminus\{e,f\})=\rho-1\) against one
representative per existing class, using \(O(nm)\subseteq O(n^2)\) bipartite
matchings. Assign column \((0,0)^\top\in\mathbb{Z}^2\) to every loop of \(M\),
and column \((1,j)^\top\in\mathbb{Z}^2\) to every element of \(C_j\)
(\(1\le j\le m\)). Any two vectors from distinct classes \(C_j\ne C_{j'}\) are
linearly independent over \(\mathbb{Q}\), while vectors in the same class are
identical, so this gives an exact \(\mathbb{Q}\)-representation of \(M\) with
\(O(\log n)\)-bit integer entries.

Applying the [represented-matroid DP](linear-source-quotas.md#dynamic-program-and-complexity)
at rank \(r=2\) gives \(m_{\text{wit}}\le 2(k+1)\),
\(N_f=\sum_{i=0}^{\min(k,m_{\text{wit}})}\binom{m_{\text{wit}}}{i}\le 2^{2k+2}\),
\(O(k)\) shared variables, and coefficient dimension
\((1+2(k+1)N_f)^{O(k)}=2^{O(k^2)}\). Since the representation entries have bit
length \(O(\log n)\), the total deterministic bit complexity is
\(2^{O(k^2)}\operatorname{poly}(n,L)\), independent of \(\rho\) and \(u\).

### Rank-\(r\) primal transversal matroids

Given a bipartite presentation of a rank-\(r\) transversal matroid \(T\) on
\(n\) ground elements, applying Gurjar et al. [5, Theorem 1] directly without
dualization produces an exact \(\mathbb{Q}\)-representation of \(T\) in
\(2^{O(\log^2(r+1))}\operatorname{poly}(n)\) time with entry bit length
\(2^{O(\log^2(r+1))}\). For \(r,k\ge 1\), running the
[represented-matroid DP](linear-source-quotas.md#dynamic-program-and-complexity)
on this representation runs in
\(2^{O(\log^2(r+1)+kr(r+k)\log(r+1))}\operatorname{poly}(n,L)=2^{O(kr(r+k)\log(r+1))}\operatorname{poly}(n,L)\)
bit operations, matching the represented exponent independently of the number
of groups and protected groups. The \(r=0\) and \(k=0\) cases run in
polynomial time as above. \(_\blacksquare\)

---

## Proof of Theorem 3: Rank-oracle lower bound

Fix \(k\ge 1\) and \(r\ge 2\), set \(t=\lfloor r/2\rfloor\ge 1\), and choose
any integer \(u\ge t+1\). Put \(M=2t+5\). Define a ground set \(E\) with \(u\)
disjoint groups \(G_1,\ldots,G_u\), each containing:
- one safe element \(s_i\) of cost \((k+1)M\),
- two vulnerable elements \(a_i,b_i\) of cost \(M+1\) each, and
- \(k-1\) vulnerable filler elements \(F_i=\{f_{i,1},\ldots,f_{i,k-1}\}\) of
  cost \(M\) each (empty for \(k=1\)).

Add an ungrouped set \(Z=\{z_1,\ldots,z_k\}\) of \(k\) vulnerable elements of
cost \(0\), and if \(r=2t+1\) is odd, add one ungrouped safe element \(c_0\) of
cost \(0\). Set \(S_{\text{core}}=\bigl(\bigcup_{i=1}^u\{a_i,b_i\}\bigr)\cup Z\),
so \(|S_{\text{core}}|=2u+k\ge 2t+k+2>2t\). For each \(I\in\binom{[u]}{t}\),
define \(H_I=\bigcup_{i\in I}\{a_i,b_i\}\subseteq S_{\text{core}}\). Notice
that distinct \(t\)-subsets \(I\ne I'\) satisfy
\(|H_I\cap H_{I'}|=2|I\cap I'|\le 2t-2\).

For any collection \(\mathcal{H}\subseteq\{H_I:I\in\binom{[u]}{t}\}\),
define \(M_{\mathcal{H}}\) on \(E\) as the direct sum of:
1. the rank-\(2t\) sparse paving matroid on \(S_{\text{core}}\) with
   circuit-hyperplane set \(\mathcal{H}\), whose rank function on
   \(Q\subseteq S_{\text{core}}\) is \(2t-1\) if \(Q\in\mathcal{H}\) and
   \(\min(|Q|,2t)\) otherwise;
2. loops on all elements \(s_i\) and \(F_i\); and
3. the rank-\(1\) free matroid on \(\{c_0\}\) when \(r\) is odd.

Because every pair of distinct sets in \(\mathcal{H}\) intersects in at most
\(2t-2\) elements, \(M_{\mathcal{H}}\) is a valid rank-\(r\) matroid for
every \(\mathcal{H}\). (For the loopless variant at \(r\ge 3\), take
\(t=\lfloor(r-1)/2\rfloor\ge 1\), place a safe element \(c_{\text{lp}}\) of
cost \(0\) together with all \(s_i\) and \(F_i\) into a rank-\(1\) parallel
class, and add \(c_0\) iff \(r-1\) is odd; every singleton has rank \(1\), the
total rank is \(r\), and adding \(c_{\text{lp}}\) and \(c_0\) at zero cost
saturates both auxiliary summands.)

In every group \(G_i\), the cheapest safe cost is \((k+1)M\) and the \(k+1\)
cheapest vulnerable elements cost \(2(M+1)+(k-1)M=(k+1)M+2\), so all \(u\)
groups are protected and any feasible purchase \(P\) spends at least \((k+1)M\)
on \(P\cap G_i\). Consider \(\mathcal{H}_{\text{all}}=\{H_I:I\in\binom{[u]}{t}\}\)
versus \(\mathcal{H}_{I^*}=\mathcal{H}_{\text{all}}\setminus\{H_{I^*}\}\) for a
hidden \(t\)-subset \(I^*\in\binom{[u]}{t}\).

- **Upper bound for \(M_{\mathcal{H}_{I^*}}\).** Select
  \(P^*=\bigcup_{i\in I^*}(\{a_i,b_i\}\cup F_i)\cup\bigcup_{i\notin I^*}\{s_i\}\cup Z\)
  (together with \(c_0\) and \(c_{\text{lp}}\) when present). Its cost is
  \((k+1)Mu+2t\). Every group satisfies its accessibility quota. Moreover
  \(|P^*\cap S_{\text{core}}|=|H_{I^*}\cup Z|=2t+k\), so deleting at most
  \(k\) vulnerable elements leaves at least \(2t\) elements in
  \(H_{I^*}\cup Z\). Every set \(H_I\in\mathcal{H}_{I^*}\) contains a pair
  \(\{a_j,b_j\}\) with \(j\notin I^*\) and therefore satisfies
  \(H_I\not\subseteq H_{I^*}\cup Z\).
  Thus every surviving set has full rank \(r\), and
  \(\mathrm{OPT}(M_{\mathcal{H}_{I^*}})\le (k+1)Mu+2t\).

- **Lower bound for \(M_{\mathcal{H}_{\text{all}}}\).** Suppose a purchase
  \(P\) has cost at most \((k+1)Mu+2t+1\). Because \(M=2t+5>2t+1\), no group
  \(G_i\) can buy \(s_i\) together with any additional element of \(G_i\). Hence
  each group \(G_i\) buys either \(\{s_i\}\) (extra cost \(0\)) or
  \(\{a_i,b_i\}\cup F_i\) (extra cost \(2\)), and at most \(t\) groups \(I\)
  can choose the vulnerable option. Failing the \(k\) vulnerable elements of
  \(Z\) leaves \((P\setminus Z)\cap S_{\text{core}}\subseteq H_I\), which has
  size \(2|I|\le 2t\) and rank at most \(2t-1\) in \(M_{\mathcal{H}_{\text{all}}}\)
  because \(H_I\in\mathcal{H}_{\text{all}}\) whenever \(|I|=t\). On the other
  hand, selecting \(\{a_i,b_i\}\cup F_i\) in \(t+1\) groups, \(\{s_i\}\) in the
  remaining groups, and \(Z\) (plus \(c_0,c_{\text{lp}}\)) leaves at least
  \(2t+2\) elements in \(S_{\text{core}}\) after any \(k\) failures, giving
  \(\mathrm{OPT}(M_{\mathcal{H}_{\text{all}}})=(k+1)Mu+2t+2\).

Because \(r_{M_{\mathcal{H}_{\text{all}}}}(Q)=r_{M_{\mathcal{H}_{I^*}}}(Q)\)
for every \(Q\subseteq E\) with \(Q\cap S_{\text{core}}\ne H_{I^*}\), any
deterministic algorithm making fewer than \(\binom{u}{t}\) rank or
independence queries returns the same answer on \(M_{\mathcal{H}_{\text{all}}}\)
and some \(M_{\mathcal{H}_{I^*}}\), contradicting the optimum separation.
For randomized algorithms succeeding with probability at least \(2/3\) using
\(q\) queries, draw \(I^*\in\binom{[u]}{t}\) uniformly; an input from
\(\mathcal{H}_{I^*}\) follows the same execution transcript as
\(M_{\mathcal{H}_{\text{all}}}\) unless a query hits \(Q\cap S_{\text{core}}=H_{I^*}\)
(probability at most \(q/\binom{u}{t}\)), so Yao's minimax bound gives
\(\frac{2}{3}\le\frac{1}{3}+\frac{q}{\binom{u}{t}}\), hence
\(q\ge\frac{1}{3}\binom{u}{t}\). \(_\blacksquare\)

---

## Proof of Theorem 4 and Proposition 5: Strict Gammoid \(q\)-Parity and tensor obstruction

### Reduction from Strict Gammoid \(q\)-Parity (Theorem 4)

Let \((D=(V,E_D),T,G_1,\ldots,G_g)\) be an instance of Strict Gammoid
\(q\)-Parity with \(q\ge 3\), \(|T|=d=qt\) (\(t\ge 1\)), and pairwise
disjoint \(q\)-element groups \(G_i=(v_{i,1},\ldots,v_{i,q})\subseteq V\).
Strip any self-loops from \(E_D\), as self-loops play no role in vertex-disjoint
paths. Form the augmented digraph \(D_1=(V_1,E_1)\) by adding a disjoint sink
copy \(T'=\{v':v\in T\}\) and \(q-1\) vertices \(Z=\{z_1,\ldots,z_{q-1}\}\) to
\(V\), with
\[
  E_1 = E_D \cup \{(v,v'):v\in T\} \cup \{(z_j,v'):z_j\in Z,\ v'\in T'\}.
\]
Because every \(v\in T\) has unique out-neighbor \(v'\in T'\), a subset
\(X\subseteq V\) is linkable to \(T'\) in \(D_1\) iff it is linkable to
\(T\) in \(D\). Moreover, for any \(X\subseteq V\) and \(F_Z\subseteq Z\),
directing \(Z\setminus F_Z\) to unused sinks of \(T'\) gives
\[
  r_{M(D_1,T')}(X\cup(Z\setminus F_Z)) = \min\bigl(d,\,r_{M(D,T)}(X)+|Z\setminus F_Z|\bigr).
\]

Now construct the directed state graph \(G^*=(V^*,E^*)\) as follows:
- **Physical and center states:** Include one state for each \(v\in V_1\) and
  one center state \(c_v\) for each \(v\in V_1\setminus T'\). Add arcs
  \((v,c_v)\), \((c_v,v)\), and \((c_v,u)\) for every \((v,u)\in E_1\).
  Vertices in \(T'\) have outdegree zero.
- **Source group cycles:** For each \(i\in\{1,\ldots,g\}\), introduce one safe
  state \(s_i\), \(q\) vulnerable copy states \(a_{i,1},\ldots,a_{i,q}\), and
  \(q+1\) vulnerable filler states \(f_{i,0},f_{i,1},\ldots,f_{i,q}\) arranged
  in the directed cycle
  \[
    s_i\to f_{i,0}\to a_{i,1}\to f_{i,1}\to a_{i,2}\to f_{i,2}\to\cdots\to a_{i,q}\to f_{i,q}\to s_i,
  \]
  together with port arcs \((f_{i,r-1},v_{i,r})\) for \(1\le r\le q\).
- **Loop-only connector:** Add arcs \((f_{1,q},c_v)\) for all
  \(v\in V_1\setminus T'\).

All vertices outside \(\{s_1,\ldots,s_g\}\) are vulnerable. Set \(k=q-1\),
\(M_{\text{big}}=2d+1\), \(\tau=qM_{\text{big}}g+d\), and weights:
\[
  w(s_i)=qM_{\text{big}},\quad
  w(a_{i,r})=M_{\text{big}}+1,\quad
  w(z_j)=0,\quad
  w(x)=2(qM_{\text{big}}+d+1)\ \text{for all remaining }x.
\]

**Structure of \(G^*\).** The source SCCs of \(G^*\) are precisely the \(g\)
group cycles \(C_i=\{s_i\}\cup\{a_{i,r}\}\cup\{f_{i,r}\}\): their vertices
receive arcs only from within their own cycles, whereas every \(c_v\) receives
from \(f_{1,q}\), every \(v\in V_1\setminus T'\) receives from \(c_v\), and
every \(v'\in T'\) receives from \(c_{z_1}\) (as \(q-1\ge 2\)).

In the bipartite state graph \(B_{G^*}\), the row set
\(R_{\text{loop}}=(V_1\setminus T')\cup\{s_i\}\cup\{a_{i,r}\}\cup\{f_{i,q}\}\)
and column set
\(C_{\text{loop}}=\{c_v:v\in V_1\setminus T'\}\cup\{f_{i,r}\}\cup\{s_i\}\)
form a disconnected square block (no edges cross between
\((R_{\text{loop}},C_{\text{loop}})\) and its complement) saturated by the
perfect matching \(v\leftrightarrow c_v\), \(s_i\leftrightarrow f_{i,0}\),
\(a_{i,r}\leftrightarrow f_{i,r}\) (\(1\le r\le q\)), and
\(f_{i,q}\leftrightarrow s_i\). Hence every vertex in \(C_{\text{loop}}\) is a
coloop of \(T(B_{G^*})\) and a loop of \(M^*=T(B_{G^*})^*\). Excluding the
\(d\) isolated rows \(T'\), the complementary bipartite block
\(B_{\text{rem}}\) on rows \(\{c_v:v\in V_1\setminus T'\}\cup\{f_{i,r-1}\}\)
and columns \(V_1\cup\{a_{i,r}\}\) satisfies \(T(B_{G^*})^*|_{V^*\setminus C_{\text{loop}}}=T(B_{\text{rem}})^*\);
in \(B_{\text{rem}}\), matching \(c_v\leftrightarrow v\) and
\(f_{i,r-1}\leftrightarrow a_{i,r}\) leaves \(T'\) unmatched, so \(G^*\) has
matching deficiency \(|T'|=d\). By Lemma 1, the associated strict-gammoid
digraph coincides with \(D_1\) on \(V_1\) with an arc \(a_{i,r}\to v_{i,r}\)
from each \(a_{i,r}\). Thus every \(a_{i,r}\) is a parallel copy of \(v_{i,r}\)
in \(M^*\) and \(M^*|_{V_1}=M(D_1,T')\).

**Equivalence of decisions.**
- Suppose \(J\subseteq [g]\) with \(|J|=t\) makes \(\bigcup_{i\in J}G_i\) a
  basis of \(M(D,T)\). Select
  \(P=\{a_{i,r}:i\in J,\,1\le r\le q\}\cup\{s_i:i\notin J\}\cup Z\). Its cost
  is \(w(P)=t\cdot q(M_{\text{big}}+1)+(g-t)qM_{\text{big}}=qM_{\text{big}}g+d=\tau\).
  Every group \(C_i\) retains at least one survivor under \(k=q-1\) failures.
  If \(F\subseteq P\cap V_{\text{vuln}}\) deletes \(m_1\) elements of
  \(\{a_{i,r}:i\in J\}\) and \(m_2\) elements of \(Z\) with
  \(m_1+m_2\le q-1\), then \(|Z\setminus F|=q-1-m_2\ge m_1\) and the
  surviving elements of \(\{a_{i,r}:i\in J\}\) project to a subset of size
  \(d-m_1\) of the basis \(\bigcup_{i\in J}G_i\) of \(M(D,T)\). By the rank
  formula for \(M(D_1,T')\), the survivors span rank
  \(\min(d,(d-m_1)+|Z\setminus F|)=d\) in \(M^*\).
- Conversely, let \(P\subseteq V^*\) be feasible of cost \(w(P)\le\tau\).
  Every source cycle \(C_i\) must survive \(k=q-1\) failures, so \(P\cap C_i\)
  contains \(s_i\) (cost \(qM_{\text{big}}\)) or at least \(q\) vulnerable
  vertices of \(C_i\) (cost at least \(q(M_{\text{big}}+1)\)). In either case
  \(w(P\cap C_i)\ge qM_{\text{big}}\), and if \(P\cap C_i\) contains any
  \(f_{i,r}\) then
  \(w(P\cap C_i)\ge 2(qM_{\text{big}}+d+1)>qM_{\text{big}}+d\). Summing over
  \(i=1,\ldots,g\) prevents \(P\) from containing any \(f_{i,r}\), any
  \(c_v\), or any \(v\in V_1\setminus Z\), as each such element has weight
  \(2(qM_{\text{big}}+d+1)\) outside the baseline sum \(qM_{\text{big}}g\).
  Thus \(P\subseteq\bigcup_{i=1}^g(\{s_i\}\cup\{a_{i,1},\ldots,a_{i,q}\})\cup Z\).
  Let \(J=\{i\in[g]:\{a_{i,1},\ldots,a_{i,q}\}\subseteq P\}\). Each \(i\in J\)
  pays excess \(q\), and any stray \(a_{i,r}\) in \(i\notin J\) adds
  \(M_{\text{big}}+1>d\); hence \(q|J|\le d=qt\) (so \(|J|\le t\)) and no
  stray \(a_{i,r}\) is bought. Applying the legal \(k\)-failure \(F=Z\) leaves
  only the loops \(\{s_i:i\notin J\}\) and \(\{a_{i,r}:i\in J,\,1\le r\le q\}\),
  forcing \(r_{M(D,T)}(\bigcup_{i\in J}G_i)=d\), hence \(|J|=t\) and
  \(\bigcup_{i\in J}G_i\) is a basis of \(M(D,T)\).

Strict Gammoid \(3\)-Parity is NP-complete by a direct reduction from
\(3\)-Dimensional Matching, admits a randomized \(\operatorname{FPT}(d)\)
algorithm via random rational representation and deterministic linear
representative sets [3], and admits a deterministic
\(2^{O(\log^2(\rho+1)+d)}\operatorname{poly}(n)\) algorithm via Gurjar et al.
[5], whereas deterministic \(\operatorname{FPT}(d)\) remains open. \(_\blacksquare\)

### Multilinear tensor obstruction (Proposition 5)

*Proof of Proposition 5.* Let \(e_1,\ldots,e_D\) be any basis of \(V\). Every
elementary tensor \(e_{i_1}\otimes\cdots\otimes e_{i_s}\in V^{\otimes s}\) and
every elementary exterior product \(e_{i_1}\wedge\cdots\wedge e_{i_s}\in\bigwedge^s V\)
uses at most \(s\) basis vectors, so \(W=\operatorname{span}\{e_{i_1},\ldots,e_{i_s}\}\)
satisfies \(\dim W\le s\le (d-1)t\) and contains that basis element in
\(W^{\otimes s}\) (respectively \(\bigwedge^s W\)). Since elementary tensors
span \(V^{\otimes s}\) and elementary exterior products span \(\bigwedge^s V\),
both spans equal the full ambient space. \(_\blacksquare\)

In Koana and Kumabe [8, Theorem 1.5], representative families of size \(2r^d\)
are computed in \(N^{O(d)}\) time for families of \(d\)-tuples across a
direct sum of \(d\) matroids, because each coordinate \(i\in[d]\) has its own
quotient map \(\pi_i:\mathcal{V}_i\to\mathcal{V}_i/S_i(Y_i)\) and bad tuples
map into \(\sum_{i=1}^d\pi_i(S_i(C_i))\otimes\bigotimes_{j\ne i}(\mathcal{V}_j/S_j(Y_j))\).
Compressing partial subsets of size \(s<d\) inside a single scale-\(t\)
representation of rank \(d\) via subspaces of \(V^{\otimes s}\) or
\(\bigwedge^s V\) encounters Proposition 5: the family of
all rank-deficient subspaces \(W\) of dimension at most \((d-1)t\) spans the
entire ambient space, whereas passing to \(\bigwedge^{st}V\) requires dimension
\(\binom{dt}{st}=2^{\Omega(n^2)}\) since \(t\ge 2n^2/\delta\).

---

## Reproduction and references

The deterministic algorithms, cut-covering verifier, sparse-paving rank-oracle
construction, and Strict Gammoid \(q\)-Parity reduction are implemented in
[scripts/deterministic_actuator_placement.py](scripts/deterministic_actuator_placement.py).
Run `python scripts/audit_deterministic_actuator_placement.py` from the repository
root (Python 3.10+, standard library only):

- **275** exhaustive `d = 0, 1, 2` state-graph optimum comparisons (95 at `d = 0`,
  95 at `d = 1`, 85 at `d = 2`, verifying exact rank-two integer representations
  on all \(2^n\) subsets) plus **15** larger `d = 1` graphs (`n = 18`, `u = 5`,
  `k = 3`);
- **65** full \(2^n\)-subset Ingleton–Piff rank-table checks against an
  independent vertex-disjoint max-flow oracle, **30** mimicking digraph
  constructions from exact cut covers (**27** strictly reducing \(\rho_Z<\rho\))
  verifying \(T(B_Z)^*|_S=M|_W\) and compressed-optimum preservation, and **2**
  exact finite-field \(\mathbb{F}_p\) checks of Koana–Kumabe Wronskian matrices
  and Lemma 4.5 dualization;
- **30** sparse-paving oracle instances (**299,776** exhaustive matroid rank
  axiom checks and **22** indistinguishable-pair exact optimum separations across
  loops and loopless variants);
- **14** Strict Gammoid \(3\)-Parity actuator reductions verifying matching
  deficiency \(d\), exact source SCCs, parallel/loop structure, restricted dual
  matroid equality against \(M(D_1,T')\), and threshold feasibility equivalence.

1. A. W. Ingleton and M. J. Piff, [*Gammoids and transversal matroids*](https://doi.org/10.1016/0095-8956(73)90031-2), *Journal of Combinatorial Theory, Series B* 15 (1973), 51–68.
2. E. L. Lawler, [*Matroid intersection algorithms*](https://doi.org/10.1007/BF01681329), *Mathematical Programming* 9 (1975), 31–56.
3. F. V. Fomin, D. Lokshtanov, F. Panolan, and S. Saurabh, [*Efficient computation of representative families with applications in parameterized and exact algorithms*](https://doi.org/10.1145/2886094), *Journal of the ACM* 63(4) (2016), 29:1–29:60.
4. S. Pequito, S. Kar, and A. P. Aguiar, [*A framework for structural input/output and control configuration selection in large-scale systems*](https://doi.org/10.1109/TAC.2015.2437525), *IEEE Transactions on Automatic Control* 61 (2016), 303–318.
5. R. Gurjar, D. Lokshtanov, P. Misra, F. Panolan, S. Saurabh, and M. Zehavi, [*Quasipolynomial-time deterministic kernelization and (gammoid) representation*](https://doi.org/10.4230/LIPIcs.MFCS.2025.54), *MFCS 2025*, LIPIcs 345, 54:1–54:17.
6. A. Chatterjee, S. Ghosh, R. Gurjar, R. Raj, and T. Thierauf, [*Bipartite matching is in NC*](https://eccc.weizmann.ac.il/report/2026/100/), ECCC TR26-100, 2026.
7. S. Kopparty and S. Saraf, [*On the CGGRT criterion for detecting bipartite perfect matchings in NC*](https://arxiv.org/abs/2607.15554), arXiv:2607.15554v1, 2026.
8. T. Koana and S. Kumabe, [*A deterministic polynomial kernel for Odd Cycle Transversal*](https://arxiv.org/abs/2609.29141), arXiv:2609.29141v1, 2026.
9. M. Kundu, S. Nandi, and S. Saurabh, [*When elements fail: Weighted fault-tolerant matroid bases*](https://arxiv.org/abs/2611.16333), arXiv:2611.16333v1, 2026.
