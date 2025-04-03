## name
mvTreeScale
## title
An Markov chain Monte Carlo (MCMC) move for Time tree scaling 
## description
mvTreeScale is a Markov chain Monte Carlo (MCMC) move that scales the entire TimeTree through the root age while keeping the tree toplogy intact
## details
The mvTReeScale move applies to a given TimeTree, adjusting all node ages while modifying the root age 
## authors
## see_also
mvSubtreeScale
mvNodeTimeSlide
mvScale
## example

#define a simple TimeTree
tree ~ dnBDP(lambda=1.0, mu=0.2, rootAge=10, samplingFraction=0.8)
#Assign a TreeScale move to the tree
moves.append(mvTreeScale(tree=tree, delta=1, tune=TRUE, weight=5))


## references
