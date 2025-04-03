## name
mvSPR
## title
Subtree Pruning and Regrafting (SPR) Move for Tree Rearrangement
## description
mvSPR is a tree topology move in Markov Chain Monte Carlo (MCMC) sampling
## details
mvSPR is an MCMC move that changes tree topology using Subtree Pruning and Regrafting (SPR). It cuts a subtree and reattaches it elsewhere, allowing larger jumps in tree space than Nearest-Neighbor Interchange (NNI), helping MCMC explore different topologies more efficiently
## authors
## see_also
mvNNI
## example
topology ~ dnUniformTopology(taxa, outgroup=out_group)
moves.append( mvSPR(topology, weight=ntaxa*2.0) )
## references
