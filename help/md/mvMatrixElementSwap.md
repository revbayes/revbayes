## name
mvMatrixElementSwap
## title
Matrix Element Swap Move
## description
Proposes exchanging two elements of a matrix.
## details
The `mvMatrixElementSwap` move picks two elements of a matrix at random and
exchanges them. The proposal only rearranges the value, so it is symmetric. It
is useful where the elements are exchangeable and a move on one element at a
time mixes poorly between their orderings.

By default any two elements may be exchanged. The swap can instead be confined
to one row or one column, which keeps it meaningful when the rows or columns of
a matrix hold different quantities. The `margin` argument follows R's MARGIN:
1 swaps within a randomly chosen row, 2 within a randomly chosen column. The
`row` and `col` arguments name the row or column instead, counting from 1.
Only one of `margin`, `row` and `col` may be given.
## authors
June Walker
## see_also
mvMatrixElementScale
mvMatrixElementSlide
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any two elements of the matrix
    moves.append( mvMatrixElementSwap(x, weight=1) )
    # two elements of column 2
    moves.append( mvMatrixElementSwap(x, col=2, weight=1) )

## references
