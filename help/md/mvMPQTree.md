## name
mvMPQTree

## title
Tree proposals for reversible and non-reversible substitution analyses

## description
Proposes tree updates in analyses using the reversible and non-reversible
nucleotide models implemented by `dnQ`, optionally updating root position.

## details
`tree` should be a stochastic, non-clock tree parameterized by branch
lengths, as produced by `dnUniformTopologyBranchLength`.

`updateRoot` defaults to TRUE and allows proposals that change root
position. Set it to FALSE when the root must remain on a specified
outgroup branch. This move does not automatically detect an assigned
outgroup or enforce that root constraint.

Root position can affect the likelihood under a non-reversible model.
Under the stationary time-reversible model, the likelihood is unchanged
by moving the root on a fixed unrooted tree with fixed branch lengths;
the root must nevertheless be handled consistently with the tree prior.

This move does not take a rate-matrix argument. Connect the tree and
the rate matrix through the phylogenetic likelihood, using `dnPhyloCTMC`.

## authors
John P. Huelsenbeck
Sebastian Höhna

## see_also
dnQ
dnUniformTopologyBranchLength
dnPhyloCTMC
mvMPQRateMatrix
mvReversibilityAwareBranchLengthScale

## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"))

    tree ~ dnUniformTopologyBranchLength(
        taxa=taxa, branchLengthDistribution=dnExp(10.0), rooted=TRUE)

    moves = VectorMoves()
    moves.append(mvMPQTree(tree=tree, updateRoot=TRUE, weight=4))

    # If the root must stay on an assigned outgroup branch, instead
    # construct the move with updateRoot=FALSE.
