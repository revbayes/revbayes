## name
mvScale
## title
Proportional Scaling Move for MCMC Simulations
## description
The mvScale move is proposal mechanism in Markov Chain Monte Carlo (MCMC) simulation. It proposes multiplicative updates to continuous parameters constrained to be positive, such as rates or branch lengths in phylogenetic models.
## details
The mvScale move updates a parameter by multiplying it with a randomly chosen factor from a proposal distribution. The size of this factor is determined by lambda, which affects how likely proposed changes are to be accepted 
## authors
## see_also
mvSlide
mvTreeScale
## example
speciation_rate ~ dnExponential(10)
moves.append( mvScale(speciation_rate, weight=1) )
## references
