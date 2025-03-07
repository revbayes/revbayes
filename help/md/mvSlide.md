## name
mvSlide
## title
Sliding Window Move for Continuous Parameters
## description
The mvSlide function is sliding-window moves propose an update by drawing a random number from a uniform distribution and then adding this random number to the current value
## details
The move have two arguments weight and delta. The argument weight decides how often this move is used in MCMC. The argument delta controls the size of proposed changes.
## authors
## see_also
mvScale 
move_uniform
## example
p ~ dnUniform(0,1)g
moves.append(mvSlide(p, delta=0.05, weight=1))
## references
