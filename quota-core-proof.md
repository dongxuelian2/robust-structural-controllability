# Unit-cost robust spanning with arbitrarily many disjoint quotas

## Theorem

Let M be a rank-d matroid on E, supplied by an independence oracle. Let S_1,...,S_s be pairwise disjoint subsets, and put G_0=E minus their union. All selected elements may fail. For k>=0, a minimum-cardinality set D satisfying

- rank(D minus F)=d for every F subset D with |F|<=k, and
- |D intersect S_i|>=k+1 for every i

can be found in f(d,k) times a polynomial in the input and oracle cost. The parameter does not include s. The resulting bounded set preserves a robust core, not necessarily the whole quota-completed solution.

Check first that |S_i|>=k+1 for every i; otherwise the instance is infeasible. When d=0, taking any k+1 elements of each S_i is optimal. Henceforth d>=1. Set b=k+1, q=bd, and

    H=binom(d+k-1,k).

This follows from the Bollobas set-pairs inequality. For each deficient
hyperplane P_i=cl(Y-F_i), choose a (d-1)-element basis A_i in Y intersect P_i,
and let D_i=Y minus P_i, where |D_i|<=k. Distinct hyperplanes yield
cross-intersecting pairs (A_i,D_i), giving at most binom(d+k-1,k) choices.

## Robust cores and the completion objective

Every robust spanning set D contains a robust spanning subset C of size at most q. To see this, successively choose k+1 disjoint maximal independent sets of the remaining restriction of D, and take their union C. Each omitted element belongs to the closure of each chosen layer. For any at most k failures in C, one of the layers avoids all failures, so D minus F is contained in cl(C minus F). Robustness of D implies robustness of C. Loops can be omitted from C.

For a robust core C define

    phi(C)=|C intersect G_0|+sum_i max(0,|C intersect S_i|-b).

The cheapest quota completion of C has cardinality exactly bs+phi(C): for each S_i with fewer than b selected elements, add any missing number of elements from S_i. Robust spanning is monotone under adding elements: for any failure set in the completed set, the surviving core corresponds to at most k failures in C and therefore spans.

Consequently the original optimum is bs plus the minimum phi(C) over robust cores of size at most q. Indeed, a small core inside any feasible D has phi(C)<=phi(D)=|D|-bs, while completion of any robust core is feasible. Original ground-set elements, including loops, remain available for completion.

An ordinary robust-spanning core need not preserve this objective. In the rank-one uniform matroid on {a,b,c,d}, let the only quota-two group be {c,d}. A core retaining only {a,b} preserves the minimum robust-spanning cardinality two, but its best quota completion has size four; the original quota optimum {c,d} has size two. The color argument below prevents this loss.

## Exchange lemma

Fix A subset E containing no loops, with rank(A)=h>=1. Let Z subset A be h-uniform, meaning every subset of size at most h is independent. Suppose |Z|>=h, so Z spans A. Let B be any robust spanning set of size at most q and x in B intersect A. Write Y=B minus {x}.

For each F subset Y with |F|<=k for which rank(Y minus F)<d, robustness of B implies rank(Y minus F)=d-1 and x lies outside the hyperplane cl(Y minus F). There are at most H distinct such hyperplanes: each has a (d-1)-element basis contained in Y, whose size is at most q-1. Every such hyperplane contains at most h-1 elements of Z, because any h elements of Z span A and hence would span x.

Therefore at most (h-1)H elements of Z are forbidden by deficient failure hyperplanes, and at most q-1 more belong to Y. Any remaining y yields a robust set B' = Y union {y}. For failures not containing y, y repairs every deficient hyperplane. For failures containing y, replace y in the failure set by x; the surviving sets of B and B' coincide. In all applications below x is outside Z, so these bounds also ensure y is outside B.

Define

    c_h=(h-1)H+q,       T_h=(q+1)c_h.

Any subset of Z of size c_h therefore contains a valid exchange candidate.

## One capped-uniform recursion

Treat G_0,S_1,...,S_s as colors, including empty G_0 if desired. Define Core(A), where A contains no loops and h=rank(A):

1. Return empty if A is empty. Returning all of A when |A|<=T_h is optional.
2. Scan A in any fixed order, maintaining Z initially empty. Accept x exactly when adding x preserves h-uniformity and its color currently has fewer than c_h accepted elements. Stop once |Z|=T_h or after exhausting A.
3. If h=1, return Z. Otherwise, for every (h-1)-subset T of Z, recursively compute Core(A intersect cl(T)). Return Z together with all these child cores.

Step 3 is required in both stopping cases. No weight ordering is needed.

The final Z contains at least h elements. Before reaching h elements, the cap cannot bind because c_h>=q>=d>=h, and uniformity testing is ordinary greedy independence. Thus exhaustion before h would contradict rank(A)=h. Each child has rank exactly h-1 because it contains the independent set T and lies in cl(T).

## Strong preservation property

For every robust spanning B of size at most q, there is robust spanning B' of the same size with

    B' minus A = B minus A,
    B' intersect A subset Core(A),
    phi(B') <= phi(B).

Prove this by induction on h. Define U=Z when h=1, and otherwise let U be the union of the child ground sets A intersect cl(T), not the union of their output cores. For h>=2, Z is contained in U.

Every rank-rejected scanned element belongs to U. A dependence witnessing rejection contains that element and at most h-1 elements of Z, and those retained elements span it. Pad them to an (h-1)-subset of final Z. For h=1 no nonloop element can be rank-rejected.

Take x in (B intersect A) minus U. There are two cases.

**Scanned x.** It must have been rejected by the color cap. Its color has exactly c_h representatives in final Z. The exchange lemma supplies a robust exchange into one of them. This preserves all color counts, and therefore preserves phi exactly.

**Unscanned x.** The scan reached T_h. Set Y=B minus {x}. Forbid G_0 and every actual source color S_i with |Y intersect S_i|>=b. There are at most q forbidden colors: at most q-1 actual saturated colors, plus G_0. Each has at most c_h representatives in Z. At least T_h-q c_h=c_h representatives remain. The exchange lemma supplies an unused valid y among them. Its source color has fewer than b elements in Y, so adding y has zero marginal phi. Deleting x never increases phi. Hence this exchange does not increase phi.

Each exchange inserts into Z subset U and removes an element outside U. After finitely many exchanges, the A-part lies in U; size and robustness are maintained, as is the outside-A part. This proves the h=1 case directly.

For h>=2, process each child using its induction hypothesis. Let W=Core(A) and Bad(B)=(B intersect A) minus W. A child operation changes elements only inside that child and replaces its entire intersection by elements of the child core, which is contained in W. Thus it cannot create a new bad element, and removes every bad element in that child. Since all bad elements initially lie in the union U of the children, processing all children leaves none. Overlap among children causes no problem. Each operation preserves size, robustness, and the outside-A part and does not increase phi. The strong property follows.

## Bounds and exact algorithm

Use q and H from the original rank d throughout the recursion. Let

    L_0=0,
    L_1=T_1,
    L_h=T_h+binom(T_h,h-1)L_(h-1)  for 2<=h<=d.

Then |Core(E minus loops)|<=L_d, independently of the number of source groups. The optional small-input shortcut also respects this bound.

For an explicit coarse bound, let T=T_d. Each recursion node has at most T^(d-1) children and depth at most d. Thus there are at most N=d*T^((d-1)^2) nodes (also valid for d=1), and the output has at most NT elements. At a node, compute its rank by greedy oracle queries. For each scanned element, check the color counter and independence of the subsets containing that element and up to h-1 retained elements. There are at most sum_(j=0)^(h-1) binom(T,j) such tests. For every child, membership in its closure is tested by an independence query on T union {x}. All these sets have size at most d; initial greedy rank computation may additionally query sets of size d+1. Thus O(N*n*d*(T+1)^(d-1)) oracle calls suffice, with polynomial bookkeeping, where n=|E|. No n^d enumeration occurs.

Enumerate every subset C of the output with |C|<=q. Test robustness by enumerating at most 2^q failures and greedily computing each surviving rank. Among feasible candidates choose minimum phi(C), then complete its source quotas using original ground-set elements. There are at most (NT+1)^q candidates. If there are none, preservation proves infeasibility. All costs and bounds depend only on d,k times a polynomial in the input size.

## Application to dedicated actuators

For a state digraph let T be the transversal matroid of state vertices matchable using state edges, and M=T* its dual, of rank equal to the maximum-matching deficiency d. With dedicated selected actuators, matching feasibility is equivalent to spanning M. Accessibility after at most k actuator failures is equivalent to selecting at least k+1 vertices in every source SCC. The source SCCs are disjoint. Hence, under this standard dedicated-input model with actuator failures only and unit purchase costs, the theorem gives an exact algorithm parameterized by d+k on arbitrary directed state graphs. No strong-connectivity or bounded-source-count assumption is required.
