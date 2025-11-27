## name
mvUpDownSlide
## title
Up-Down Proposal for Joint Sliding-Window Adjustments to Multiple Parameters
## description
Simultaneously applies a sliding adjustment to multiple parameters, potentially
increasing some and decreasing others by the same amount.
## details
This move adds a random value to a set of parameters, and optionally subtracts
the same value from another set of parameters. This may improve the mixing
of parameters that we expect to be either positively or negatively correlated,
such as speciation and extinction rates or the clock rate and branch durations
in time tree inference. The parameters to be incremented and decremented can be
selected using the `.addVariable()` and `.removeVariable()` methods. The actual
value to be added or subtracted is equal to delta * (u - 0.5), where delta is
a tuning parameter and u is a random draw from the uniform distribution on
[0, 1].
## authors
## see_also
mvSlide
mvSlideBactrian
## example
    moves = VectorMoves()
    
    log_speciation_rate ~ dnNormal(-1, 0.5)
    log_extinction_rate ~ dnNormal(-1, 0.5)
    
    # Define the basic properties of the move
    delta_up_down_move = mvUpDownSlide(delta=0.05, weight=5.0)
    
    # Add variables to the move to account for their positive correlation
    delta_up_down_move.addVariable(log_speciation_rate, up=TRUE)
    delta_up_down_move.addVariable(log_extinction_rate, up=TRUE)
    
    # Apply the move
    moves.append( delta_up_down_move )

## references
