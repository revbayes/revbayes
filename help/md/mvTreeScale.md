## name
mvTreeScale
## title
An Markov chain Monte Carlo (MCMC) move for Time tree scaling 
## description
mvTreeScale scales all the ages in the TimeTree by the same amount while leaving the topology unchanged
## details
The mvTreeScale move applies to a given TimeTree, adjusting all node ages while modifying the root age 
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
- citation: Yang, Ziheng, Molecular Evolution: A Statistical Approach (Oxford, 2014; online edn, Oxford Academic, 21 Aug. 2014), https://doi.org/10.1093/acprof:oso/9780199602605.001.0001, accessed 16 Apr. 2025.
