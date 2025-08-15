## name
mvVectorScale
## title
Scaling Move on Vector of Continuous Parameters
## description
mvVectorSlide proposes additive changes to all elements of a rate vector at once
## details
A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5) and scales the current vale by a scaling factor
sf = exp( lambda * u )
where lambda is the tuning parameter of the proposal.
## authors
## see_also
mvScale
mvVectorScale 
## example
moves.append( mvVectorScale(rates, lambda=1.0) )
## references
