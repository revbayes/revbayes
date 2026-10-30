## name
mvReversibilityAwareBranchLengthScale

## title
Branch-length scaling that keeps the root edge together

## description
Scales the length of a randomly chosen non-root branch. When that branch
descends from a bifurcating root, its sibling is scaled by the same
factor, so the two root branches change together and keep their length
ratio.

## details
Use this move on a stochastic, non-clock tree parameterized by branch
lengths, such as a rooted tree drawn from
`dnUniformTopologyBranchLength`.

Every branch not attached to the root is scaled on its own, as in
`mvBranchLengthScale`. Under a stationary time-reversible substitution
model, the two root branches correspond to one edge of the unrooted
tree. Their combined length affects the likelihood, whereas their
relative lengths are informed by the prior alone. Scaling them jointly
updates the combined length without changing their ratio. The move is
also usable under the non-reversible model and does not require a
rate-matrix argument.

Use this move alongside an ordinary single-branch scaler such as
`mvBranchLengthScale`. Joint scaling alone cannot explore the ratio
of the two root-branch lengths.

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
