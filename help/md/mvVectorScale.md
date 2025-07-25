## name
mvVectorScale
## title
Scaling Move on Vector of Continuous Parameters
## description
mvVectorScale proposes additive changes to all elements of a rate vector at once
## details
When proposing a new parameters, it is important to decide how big the step size should be, too large or too small both is unaffective. A very small step, the chain will move slowly and a very large step; most proposals will get rejected. That's when scaling factore is used and to control it a tuning parameter lambda is intoduced
A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5) and scales the current vale by a scaling factor
sf = exp( lambda * u ). 
## authors
## see_also
mvScale
mvVectorScal
mvTreeScale 
## example
speciation_rate ~ dnExponential(10)
rates ~ dnDirichlet( rep(1.0, 3) )
moves.append( mvVectorScale(rates, lambda=1.0, weight=2.0) )
## references
