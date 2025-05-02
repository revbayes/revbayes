## name
mvNNI
## title
Nearest Neighbor Interchange (NNI) move.
## description
Tree topology move that performs a Nearest Neighbor Interchange (NNI) on
a rooted or unrooted tree.
## details
`mvNNI` changes tree topology by interchanging the positions of two subtrees
around an internal branch. For each selected internal edge, the move considers
two alternative topologies. As there are (n - 3) internal edges in an unrooted
tree of n taxa, every such tree has (2n - 6) NNI "neighbors" that are one NNI
move away. This neighborhood is smaller than that induced by more complex
topology moves such as `mvSPR`. As a result, `mvNNI` is computationally cheaper
than `mvSPR` and will often exhibit higher acceptance rates, but explores tree
space less thoroughly and is more likely to get stuck in local optima
(resulting in poor mixing). The RevBayes implementation of the NNI move can be
applied to both `BranchLengthTree` and `TimeTree` objects.
## authors
## see_also
mvFNPR
mvSPR
mvTreeScale
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    height ~ dnUniform(0, 10)
    moves = VectorMoves()
  
    # Apply the NNI move to an unrooted BranchLengthTree
    bltree ~ dnUniformTopology(taxa)
    moves.append( mvNNI(tree=bltree, weight=taxa.size()) )
    
    # Apply the NNI move to a rooted TimeTree
    timetree ~ dnUniformTimeTree(rootAge=height, taxa=taxa)
    moves.append( mvNNI(tree=timetree, weight=taxa.size()) )

## references
- citation: Robinson DF (1971). Comparison of labeled trees with valency three. Journal of Combinatorial Theory, Series B, 11(2):105-119.
  doi: 10.1016/0095-8956(71)90020-7
  url: https://www.sciencedirect.com/science/article/pii/0095895671900207
- citation: Waterman MS, Smith TF (1978). On the similarity of dendrograms. Journal of Theoretical Biology, 73(4):789-800.
  doi: 10.1016/0022-5193(78)90137-6
  url: https://dornsife.usc.edu/msw/wp-content/uploads/sites/236/2023/09/msw-029.pdf
