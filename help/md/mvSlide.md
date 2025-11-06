## name
mvSlide
## title
Sliding-Window Move
## description
Proposes additive updates to continuous parameters.
## details
The `mvSlide` move updates a parameter by drawing a random number from
a uniform distribution and adding the draw to the current value. The move takes
arguments that control how often it should be used (`weight`) and its "window
size", i.e., the width of the uniform distribution, which determines the size
of the changes it proposes (`delta`). When the `tune` argument is set to
`TRUE`, the value of `delta` is automatically adjusted so that the acceptance
rate of the move reaches `tuneTarget`.
## authors
## see_also
mvScale
mvSlideBactrian
## example
    moves = VectorMoves()
    p ~ dnUniform(0, 1)
    moves.append(mvSlide(p, delta=0.05, weight=1))

## references
