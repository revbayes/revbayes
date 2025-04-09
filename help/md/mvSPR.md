## name
mvSPR
## title
Subtree Pruning and Regrafting (SPR) Move for Tree Rearrangement
## description
mvSPR is a tree topology move used in Markov Chain Monte Carlo (MCMC) sampling
## details
mvSPR is an MCMC move that changes tree topology using Subtree Pruning and Regrafting (SPR). SPR cuts off a subtree and reattaches it elsewhere in the original tree using the same subtree branch that was originally cut. SPR allows larger jumps in the tree space than Nearest-Neighbor Interchange (NNI), helping MCMC explore a broader range of different topologies.
## authors
## see_also
mvNNI
mbSubtreeSwap
## example
topology ~ dnUniformTopology(taxa, outgroup=out_group)
moves.append( mvSPR(topology, weight=ntaxa*2.0) )
## references
Felsenstein, J. (2004). Inferring phylogenies (pp. 41-46). Sinauer Associates.

Wikipedia contributors. Tree rearrangement. Wikipedia, The Free Encyclopedia. Retrieved 3/26/2025, from https://en.wikipedia.org/wiki/Tree_rearrangement
