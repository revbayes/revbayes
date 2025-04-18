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
mvScaleBactrian
## example
speciation_rate ~ dnExponential(10)
extinction_rate ~ dnExponential(10)
moves.append(mvUpDownScale([speciation_rate], [extinction_rate], lambda=1.0, weight=1))
## references
- citation: Yang, Ziheng, Molecular Evolution: A Statistical Approach (Oxford, 2014; online edn, Oxford Academic, 21 Aug. 2014), https://doi.org/10.1093/acprof:oso/9780199602605.001.0001, accessed 16 Apr. 2025.
