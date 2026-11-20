# A lower bound for anchored cotransversal parity gadgets

**Theorem.** Let $T$ be a transversal matroid on $k$ pairs of label
elements and a set $A$ of fixed elements, with $k\ge3$. In the dual
matroid select all of $A$ and some labels. Suppose that:

- When exactly one label from each pair is selected, the selection is
  robust precisely for one parity of the label bits.
- A selection missing both labels of a pair is never robust, even
  when every other label is selected.

Then $|A|\ge k+2^{k-1}$. The bound is tight, and in any transversal
presentation the $2k$ labels have distinct singleton neighborhoods.
In particular three-variable parity requires seven fixed elements
in this local model.

Here robust means spanning after deletion of any one selected element.
The model conditions on selecting $A$ and imposes the two requirements
by matching alone. External quotas or selectable auxiliaries can
change the conditions.

## Independent flats and row neighborhoods

Work with the complement $X$ of the selected set. Robustness is
equivalent to $X$ being an independent flat of $T$. Write $N(Y)$ for the
union of row neighborhoods of elements in $Y$. Call a canonical
$k$-set $X$ good if it is an independent flat. Each canonical $X$ contains
one label per pair. The good sets comprise one parity because taking
complements changes all bits.

Every proper partial choice of labels extends to a good set. Consequently
every canonical set $X$ is independent. If it contained a circuit,
remove one circuit element and extend the remaining circuit labels to
a good set excluding that element. This would violate flatness.

For each pair $P_i$, the two labels are independent since a good set
containing either label excludes the other. The missing-pair condition
makes $P_i$ nonflat. Any element in its closure outside the pair must
belong to $A$. If it were a label from another pair, a good set containing
that label and one member of $P_i$ would span the other member. Here
there is no parallel pair among these three elements, since a good set
containing either of two parallel labels can omit the other. A fixed
witness parallel to a label is likewise ruled out by flatness.
Thus $P_i$ and its witness form a three-element circuit. Hall's condition
gives $|N(P_i)|=2$, with the witness supported on those same rows.

For every bad canonical set $X$, some element outside $X$ lies in its
closure. It cannot be the opposite label in any pair. Replacing that
pair's label would give a good set of size $k$ spanning the same union
and hence spanning the old label. A fixed witness cannot lie in the
closure of a proper subset of $X$, since that subset extends to a
good set. Therefore $X$ and a fixed witness form a circuit of size
$k+1$. Hall's condition gives $|N(X)|=k$.

## A good set using exactly k rows

Sum $|N(X)|$ separately over the good and the bad canonical sets.
A row contributes equally to the two sums if it occurs in both labels
of a pair. It also contributes equally if it occurs in neither label
of a pair, by flipping that bit. The only remaining rows occur in
exactly one label of every pair. Each such row contributes a difference
of at most one, and there are at most two such rows because $|N(P_i)|=2$.

The bad sum is $k2^{k-1}$. Since $2^{k-1}>2$, some good set $X$
uses exactly $k$ rows. For every pair $P_i$, its two rows cannot
both be in $N(X)$, since its fixed witness would then extend $X$
dependently. Hence the chosen label of each pair has a singleton
neighborhood, and these $k$ rows are distinct by independence.
Flipping any one bit gives a bad set using exactly $k$ rows, so the
opposite label is also a singleton on a row outside $N(X)$. Flipping
any two bits gives a good independent set. Thus the opposite labels
also have distinct rows, proving the claimed structure.

With these singleton neighborhoods a pair witness must use both rows
of its pair. A bad-set witness must use all $k$ rows of that set,
since any proper subset of its rows occurs in a good set. The $k$
pair supports and $2^{k-1}$ bad-set supports are all distinct. This
proves the lower bound. Conversely take distinct singleton labels
and one fixed element on each of those supports. Hall's condition
then gives both desired properties.

## Two fixed elements suffice when the pair condition is dropped

For $k=3$ take the rows to be the four vertices of $K_4$. The labels
are its edges, each incident to its endpoints, with opposite edges
forming pairs. Add two fixed columns on distinct singleton rows.
A canonical complement is either a spanning star or a triangle.
The four stars have the same parity. Each star is an independent flat
because any additional edge or singleton column yields a matchable
four-set. Every triangle contains a fixed column's row and is nonflat.
But a complement containing two opposite edges and one edge of a
different pair is a spanning path and is also an independent flat.
Its selection misses one pair and doubles another at the same cost.
`verify_parity_guards.py` checks both of these small presentations.
