## name
mvScaleBactrian
## title
Scaling Moves Proposal Follows Bactrian Distribution
## description
mvScaleBactrian is an MCMC move that scales a parameter while applying a Bactrian-distributed perturbation. 
## details
This encourages larger jumps and improves mixing efficiency, unlike a standard scaling move that uses a log-normal proposal, mvScaleBactrian creates a bimodal proposal, like the two humps of a Bactrian camel, encouraging larger jumps
A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 and scales the current vale by a scaling factor
 sf = exp( lambda * u )
 where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
## authors
## see_also
mvScale 
## example
moves.append( mvScaleBactrian(speciation_rate_at_present,weight=5) )
## references
