## name
mvSlide
## title
Sliding Window Move for Continuous Parameters
## description
mvSlide is an MCMC move that proposes additive updates to a continuous variable
## details
mvSlide moves propose an update by drawing a random number from a uniform distribution and then adding this random number to the current value. The move have two arguments weight and delta. The argument weight decides how often this move is used in MCMC. The argument delta controls the size of proposed changes.
## authors
## see_also
mvScale 
move_uniform
## example
p ~ dnUniform(0,1)
moves.append(mvSlide(p, delta=0.05, weight=1))
## references
