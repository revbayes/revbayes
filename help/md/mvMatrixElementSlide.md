## name
mvMatrixElementSlide
## title
Matrix Element Sliding-Window Move
## description
Proposes an additive update to one element of a matrix.
## details
The `mvMatrixElementSlide` move picks one element of a matrix at random, draws
a random number from a uniform distribution and adds the draw to the element,
as `mvSlide` does for a single parameter. The "window size", i.e., the width of
the uniform distribution, is set by `delta`. When the `tune` argument is set to
`TRUE`, the value of `delta` is automatically adjusted so that the acceptance
rate of the move reaches `tuneTarget`.

By default the element is drawn from the whole matrix. The `row` and `col`
arguments confine the move to one row or one column, counting from 1, and
giving both pins a single element. This is useful when the rows or columns of
a matrix hold different quantities, so that one proposal size does not suit
them all.
## authors
Nicolas Lartillot
June Walker
## see_also
mvMatrixElementScale
mvMatrixElementSwap
mvSlide
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any element of the matrix
    moves.append( mvMatrixElementSlide(x, delta=0.1, weight=1) )
    # only the elements of column 2
    moves.append( mvMatrixElementSlide(x, delta=0.1, col=2, weight=1) )

## references
