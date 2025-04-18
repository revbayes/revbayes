## name
mvScaleBactrian
## title
Scaling Moves Proposal Follows Bactrian Distribution
## description
mvScaleBactrian is an MCMC move that scales a parameter while applying a Bactrian-distributed perturbation. 
## details
Proposes multiplicative changes to a positive parameter using a Bactrian kernel—a bimodal distribution centered at zero that avoids small steps. This leads to larger, more diverse proposals, improving mixing efficiency and reducing autocorrelation
## authors
## see_also
mvScale 
## example
moves.append( mvScaleBactrian(speciation_rate_at_present,weight=5) )
## references
- Citation: Z. Yang, & C.E. Rodríguez, Searching for efficient Markov chain Monte Carlo proposal kernels, Proc. Natl. Acad. Sci. U.S.A. 110 (48) 
  19307-19312, https://doi.org/10.1073/pnas.1311790110 (2013).
