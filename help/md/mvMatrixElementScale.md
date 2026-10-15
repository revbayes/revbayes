## name
mvMatrixElementScale
## title
Matrix Element Scaling Move
## description
Proposes a multiplicative update to one element of a matrix.
## details
The `mvMatrixElementScale` move picks one element of a matrix at random and
multiplies it with a randomly chosen factor, as `mvScale` does for a single
parameter. The size of the scaling factor is set by `lambda`. When the `tune`
argument is set to `TRUE`, the value of `lambda` is automatically adjusted so
that the acceptance rate of the move reaches `tuneTarget`. Since multiplicative
updates do not change the sign, the move is most often applied to matrices
whose elements are constrained to be positive.

By default the element is drawn from the whole matrix. The `row` and `col`
arguments confine the move to one row or one column, counting from 1, and
giving both pins a single element. This is useful when the rows or columns of
a matrix hold different quantities, so that one proposal size does not suit
them all.
## authors
Nicolas Lartillot
June Walker
## see_also
mvMatrixElementSlide
mvMatrixElementSwap
mvScale
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any element of the matrix
    moves.append( mvMatrixElementScale(x, lambda=0.1, weight=1) )
    # only the elements of column 2
    moves.append( mvMatrixElementScale(x, lambda=0.1, col=2, weight=1) )

## references
