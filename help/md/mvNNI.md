## Name
mvNNI: Nearest Neighbor Interchange(NNI) 
## Title
Nearest Neighbor Interchange(NNI) move for unrooted phylogenetic trees proposals
## description
MvNNI move applies a Nearest Neighbor Interchange(NNI) to a given phylogenetic tree, it proposes a small topological change to the tree by swapping branches 
## details
 The tree topology and branch legnths are stochastic nodes in the phylogenetic model, and NNI is one of the several tree topology moves, it swaps to adjacent subtree at an internal node
## authors
## see_also
mvSPR
mvRootTimeSlideUniform
mvSubtreeScale
## example
moves.append( mvNNI(topology, weight=num_taxa) )
## references

