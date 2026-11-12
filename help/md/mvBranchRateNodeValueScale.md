## name
mvBranchRateNodeValueScale
## title
Proportional scaling move of node values (for branch rate models)
## description
Proposes multiplicative updates to continuous node parameters.
## details

A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.

The `mvBranchRateNodeValueScale` move updates a single value for one of the nodes of a tree by multiplying it with a randomly chosen factor from a proposal (uniform) distribution. First, the `mvBranchRateNodeValueScale` picks a random node of the tree. Then, a random uniform number u ~ unif (-0.5,0.5) is drawn and scales the current vale by a scaling factor sf = exp( lambda * u ) where `lambda` is the tuning parameter of the proposal to influence the size of the proposals. If the node was not a tip node, then the two descendant nodes are scaled by the reciprocal of the scaling factor. This will ensure that node values rates

The move takes arguments that control how
often it should be used (`weight`) and the size of the scaling factor (`lambda`). When the `tune` argument is set to `TRUE`, the value of `lambda` is automatically adjusted so that the acceptance rate of the move reaches
`tuneTarget`. Since multiplicative updates do not change the sign, `mvBranchRateNodeValueScale` is most often applied to parameters that are constrained to be positive, such as
rates or branch lengths in phylogenetic models.
## authors
## see_also
mvScaleBactrian
mvSlide
## example
    moves = VectorMoves()
    speciation_rate ~ dnExponential(10)
    moves.append( mvScale(speciation_rate, weight=1) )

## references
