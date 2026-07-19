## name
mvMatrixElementSwap

## title
Matrix element swap move

## description
Exchanges two elements of a matrix, either anywhere in it or within a single row or column.

## details
By default any two elements may be exchanged. Confining the swap to one line instead exchanges the same quantity between two positions, which keeps the proposal meaningful when the other margin holds different things. `margin` follows R's MARGIN: 1 swaps within a row, 2 within a column. Supply `row` or `col` to pin the line rather than drawing one at random. Only one of `margin`, `row` and `col` may be given.

Note that the meaningful margin depends on the matrix. For the birth/death matrix of `dnFBDRP`, whose rows are taxa, only column swaps are meaningful: swapping within a row would exchange a taxon's own birth and death, and every such proposal is rejected.

The proposal permutes the value, so it is symmetric and carries no Hastings ratio. It is useful where a random walk on individual elements mixes poorly between orderings, for example over which range of a `dnFBDRP` matrix holds the oldest birth.

## authors
June Walker

## see_also
mvMatrixElementScale
mvMatrixElementSlide

## example
bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, origin_prior=dnUniform(9,15))
# column 1 holds the birth times
moves.append( mvMatrixElementSwap(bd, col=1, weight=taxa.size()) )

## references
