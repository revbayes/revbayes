## name
mvReversibilityAwareBranchLengthScale

## title
Joint scaling of the two root branches

## description
Scales the two branches descending from a bifurcating root together,
preserving their length ratio.

## details
Use this move on a stochastic, non-clock tree parameterized by branch
lengths, such as a rooted tree drawn from
`dnUniformTopologyBranchLength`.

Under a stationary time-reversible substitution model, the two root
branches correspond to one edge of the unrooted tree. Their combined
length affects the likelihood, whereas their relative lengths are
informed by the prior alone. Joint scaling updates the combined length
without changing their ratio. The move is also usable under the
non-reversible model and does not require a rate-matrix argument.

Use this move alongside an ordinary single-branch scaler such as
`mvBranchLengthScale`. Joint scaling alone cannot explore the ratio
of the two root-branch lengths.

`delta` controls proposal strength and defaults to one.
`tune` defaults to TRUE and enables tuning when requested during burn-in.
`weight` controls selection frequency relative to other moves.

## authors
John P. Huelsenbeck
Sebastian Höhna

## see_also
dnQ
dnUniformTopologyBranchLength
mvBranchLengthScale
mvMPQTree
mvMPQRateMatrix

## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"))

    tree ~ dnUniformTopologyBranchLength(
        taxa=taxa, branchLengthDistribution=dnExp(10.0), rooted=TRUE)

    moves = VectorMoves()
    moves.append(mvBranchLengthScale(tree, tune=FALSE, weight=4))
    moves.append(mvReversibilityAwareBranchLengthScale(
        tree, delta=1.0, tune=FALSE, weight=4))
