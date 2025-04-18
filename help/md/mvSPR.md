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
- citation: Swofford, D. L., & Olsen, G. J. (1990). Phylogeny reconstruction. In D. M. Hillis & C. Moritz (Eds.), Molecular Systematics (1st ed., pp.   
  411–501). Sunderland, MA: Sinauer Associates.
- citation: Felsenstein J (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach". Journal of Molecular Evolution. 17:368–76.
  doi: https://doi.org/10.1007/BF01734359
  url: https://link.springer.com/article/10.1007/BF01734359


