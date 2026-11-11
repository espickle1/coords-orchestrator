# The alignment problem, as a coordinate problem

Strip the subject matter away. What remains is a set of strings, a table of measurements on each string, and one spatial object. Each of these has its own index. An alignment relates the string indices. It does not replace them.

## Objects

- Strings \(S_0, S_1, \ldots, S_n\) over one alphabet. \(S_0\) is the reference. The others are homologous to it: the same symbols in the same order, plus substitutions, insertions, and deletions.
- For each string \(S_i\), a measurement table \(T_i\). Row \(k\) of \(T_i\) is the measurement of the \(k\)-th symbol of \(S_i\). The table has exactly as many rows as \(S_i\) has symbols. It stores no column number and no reference position.
- A spatial object \(C\), built from some string in the family or from a relative of that string. An element of \(C\) is addressed by a triple \((\text{chain},\; a,\; \text{part})\), where \(a\) is an author index assigned when \(C\) was built, and the part is a named piece of that element. A geometric query takes one element and a radius and returns a set of elements of \(C\).
- A second measurement \(E\) on \(C\), retrieved with the same triple. \(E\) is not a row of any \(T_i\).

## Four indices that look like one integer

| Index | Counts | Belongs to |
|---|---|---|
| Reference position \(p\) | Symbols of \(S_0\), from the start of \(S_0\) | \(S_0\) only |
| String position \(q_i\) | Symbols of \(S_i\), from the start of \(S_i\) | \(T_i\) |
| Column \(c\) | Columns of a multiple alignment, including gaps | The alignment |
| Author index \(a\) | Elements of \(C\), in the order \(C\) was numbered | \(C\) and \(E\) |

These four integers are equal in special cases. They are not the same coordinate.

## What an alignment supplies

A multiple alignment places every string into columns. In each column a string contributes one symbol or a gap.

Designate \(S_0\). Counting only the columns where \(S_0\) has a symbol produces the reference positions. The alignment then compiles, for each other string, a partial function

\[
f_i : \{\text{positions of } S_0\} \rightharpoonup \{\text{positions of } S_i\}.
\]

\(f_i(p) = q\) means the symbol at position \(p\) of the reference occupies the same column as the symbol at position \(q\) of \(S_i\). \(f_i(p)\) is undefined when \(S_i\) has a gap in that column. There is then no symbol to measure, and \(T_i\) has no row for it.

Two facts follow from the definition.

1. A column number is not a reference position. Columns where \(S_0\) itself has a gap do not advance \(p\). The alignment is longer than \(S_0\) whenever the reference, or any sequence placed beside it, contains an insertion relative to the others.
2. A joint alignment of \(\{S_i\}\) does not by itself name a reference. Choosing \(S_0\) and compiling each \(f_i\) is a separate step. Until that step, the alignment is a set of columns, not a shared origin.

The function \(f_i\) is the whole content of the alignment. Applying it is a different operation from having computed it.

## What the measurement tables contain

\(T_i\) is indexed by \(q_i\), the position along \(S_i\). Reading row \(p\) of \(T_i\) returns the measurement of symbol \(p\) of \(S_i\).

That symbol is the symbol at reference position \(p\) in exactly one situation: \(f_i(p) = p\). This is true for every \(p\) on the reference itself, because \(f_0\) is the identity. It is true for another string at a particular \(p\) when that string has no gap and no extra symbol anywhere before \(p\). One insertion or deletion upstream moves every later row. The shift is different for each string, because each string has its own gaps.

So the existence of the alignment does not put \(T_i\) into the reference frame. The alignment is the dictionary. The table is still filed under the original page numbers. The read in the reference frame is

\[
T_i\bigl[f_i(p)\bigr],
\]

and it is defined only when \(f_i(p)\) is defined.

## What the spatial object contains

\(C\) is indexed by author indices. The lookup for a part is \((\text{chain},\; a,\; \text{part})\). The part name selects a piece inside an element. It selects the right element only after \(a\) has been resolved.

Passing a reference position \(p\) into that lookup treats \(p\) as an author index. Those two counts were produced by different procedures. They name the same element when the person who numbered \(C\) copied the positions of \(S_0\), and otherwise they do not. The alignment of the strings never reads \(C\). \(C\) does not store \(f_i\). Homology of the underlying strings is not a numbering.

The geometric query makes this a set, not a single index. The radius returns author indices \(a_1, \ldots, a_m\). Each one needs its own translation back to a reference position and then, through \(f_i\), to a row of each \(T_i\). A neighbor on another chain belongs to another family of strings and another alignment. The reference position of the query does not address it.

## A minimal example

\[
\begin{align*}
S_0 &: \mathtt{A\;B\;C\;D\;E} \\
S_1 &: \mathtt{A\;B\;X\;D\;E} \\
S_2 &: \mathtt{A\;-\;C\;D\;E}
\end{align*}
\]

Reference position 3 is \(\mathtt{C}\).

- \(f_1(3) = 3\). Row 3 of \(T_1\) is the measurement of \(\mathtt{X}\). Same index, different symbol, and that is the correct aligned row.
- \(f_2(3) = 2\). Row 3 of \(T_2\) is the measurement of \(\mathtt{D}\). The aligned row is row 2.
- Suppose \(C\) was numbered \(10,11,12,13,14\) over \(\mathtt{A,B,C,D,E}\). The element for reference position 3 has author index 12. A query for author index 3 asks for an element \(C\) does not have, or has under a different meaning.

The alignment is present in the example. The three reads still disagree.

## The operation that the indices have to support

One query, stated only as a reference position \(p\) and a radius:

1. Translate \(p\) to an author index \(a\) on \(C\). If that translation is undefined, \(p\) is not on \(C\).
2. Run the geometric query at \(a\), producing author indices \(a_1, \ldots, a_m\).
3. Translate each \(a_j\) back to a reference position \(p_j\), then to \(f_i(p_j)\) on every string. Where either translation is undefined, record a missing element.
4. Read \(E\) at \((\text{chain},\; a,\; \text{part})\) and read \(T_i[f_i(p_j)]\).

After those four steps, both results are views of the same reference positions. Steps 1 and 3 are the alignment problem. The tables, the alignment, and the spatial lookup can each be finished while those two steps are still unapplied.

## What is already solved, and what is not

The alignment, once a reference is designated, solves the construction of \(f_i\). Each \(T_i\) solves the measurement of one string in that string's own order. \(C\) solves the geometric query in author order, and \(E\) solves the measurement addressed the same way.

The unsolved step is the identification. Nothing rewrites \(T_i\) into reference positions, and nothing rewrites the author indices of \(C\) into reference positions. An orchestrator is the procedure that holds a single integer in the reference frame and applies \(f_i\), and the structure map, before every read. Where a map is undefined it must leave a hole. Reading the next row fills that hole with a different element.
