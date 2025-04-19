## name
mvUpDownScale
## title
Up-Down Proposal for Joint Scaling of Multiple Parameters 
## description
Simultaneously scales multiple parameters up, and potentially down, by the same
factor.
## details
This move scales a set of parameters up by a random factor, and optionally
scales another set of parameters down by the same value. This may improve
the mixing of parameters that we expect to be either positively or negatively
correlated, such as speciation and extinction rates or the clock rate and
branch durations in time tree inference. The parameters to be scaled up and
down can be selected by the `.addVariable()` and `.removeVariable()` methods.
The actual scaling factor is equal to exp( lambda * (u - 0.5) ), where lambda
is a tuning parameter and u is a random draw from the uniform distribution on
[0, 1].
## authors
## see_also
mvScale
mvScaleBactrian
## example
    moves = VectorMoves()
    
    speciation_rate ~ dnExponential(10)
    extinction_rate ~ dnExponential(10)
    
    # Define the basic properties of the move
    up_down_move = mvUpDownScale(lambda=1.0, weight=5.0)
    
    # Add variables to the move to account for their positive correlation
    up_down_move.addVariable(speciation_rate, up=TRUE)
    up_down_move.addVariable(extinction_rate, up=TRUE)
    
    # Apply the move
    moves.append( up_down_move )

## references
- citation: Rannala B, Yang Z (2003). Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics, 164(4):1645-1656.
  doi: 10.1093/genetics/164.4.1645
  url: https://www.rannala.org/reprints/2003/Rannala2003a.pdf
