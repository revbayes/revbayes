## name
mvScale
## title
Proportional Scaling Move
## description
Proposes multiplicative updates to continuous parameters.
## details
The `mvScale` move updates a parameter by multiplying it with a randomly chosen
factor from a proposal distribution. The move takes arguments that control how
often it should be used (`weight`) and the size of the scaling factor
(`lambda`). When the `tune` argument is set to `TRUE`, the value of `lambda` is
automatically adjusted so that the acceptance rate of the move reaches
`tuneTarget`. Since multiplicative updates do not change the sign, `mvScale` is
most often applied to parameters that are constrained to be positive, such as
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
