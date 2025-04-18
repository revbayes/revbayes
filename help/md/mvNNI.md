## Name
mvNNI: Nearest Neighbor Interchange (NNI) 
## Title
Nearest Neighbor Interchange(NNI) move for unrooted phylogenetic trees
## description
mvNNI is a tree topology move that performs a Nearest Neighbor Interchange (NNI) on a given time tree. 
## details
mvNNI proposes small changes to a tree topology by interchanging the positions of two subtrees around an internal branch. For each selected internal edge, the move considers two alternative topologies. This move is commonly used in MCMC to explore tree space efficiently while maintaining relatively high acceptance rates
## authors
## see_also
mvSPR
mvTreeScale
## example
moves.append( mvNNI(tree=topology, weight=ntaxa) )
## references
- citation: Robinson DF (1971). Comparison of labeled trees with valency three. Journal of Combinatorial Theory, Series B, 11(2):105-119.
  doi: 10.1016/0095-8956(71)90020-7
  url: https://www.sciencedirect.com/science/article/pii/0095895671900207
- citation: Waterman MS, Smith TF (1978). On the similarity of dendrograms. Journal of Theoretical Biology, 73(4):789-800.
  doi: 10.1016/0022-5193(78)90137-6
  url: https://dornsife.usc.edu/msw/wp-content/uploads/sites/236/2023/09/msw-029.pdf