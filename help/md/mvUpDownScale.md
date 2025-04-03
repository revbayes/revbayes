## name
mvUpDownScale
## title
Up-Down Scaling Proposal for several parameters jointly. 
## description
mvUpDownScale is an MCMC proposal that simultaneously scales multiple parameters up and down by the same factor
## details
This proposal randomly scales all a set of parameter up while the other set of parameters
is scaled down by the same value. This should improve mixing in many cases.
The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
where u ~ Uniform(0,1)
## authors
## see_also
mvScale
mvScaleBartrian
## example
speciation_rate ~ dnExponential(10)
extinction_rate ~ dnExponential(10)
moves.append(mvUpDownScale([speciation_rate], [extinction_rate], lambda=1.0, weight=1))
## references
