## name
mvBranchLengthScale
## title
Proportional Scaling Move of Single Branch Lengths
## description
Proposes multiplicative updates to a branch length.
## details
The `mvBranchLengthScale` move updates a branch length by multiplying it with a randomly chosen factor from a uniform distribution. The move takes arguments that control how often it should be used (`weight`) and the size of the scaling factor (`lambda`). When the `tune` argument is set to `TRUE`, the value of `lambda` is automatically adjusted so that the acceptance rate of the move reaches `tuneTarget`. The `mvBranchLengthScale` selects a random branch in the topology for scaling. Note that you can restrict this to only terminal branches.
## authors
Sebastian Höhna
## see_also
mvScaleBactrian
mvScale
## example
    moves = VectorMoves()
    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvBranchLengthScale(phylogeny, weight=1) )

## references
