## name
mvSubtreeSwap
## title
Subtree Swap move.
## description
Tree topology move that performs a Subtree Swap on an unrooted tree.
## details
`mvSubtreeSwap` randomly picks two non-nested clades and exchanges/swaps them.
In that sense, it is a superset of the NNI move but different from an SPR move.
## authors
Sebastian Höhna
## see_also
mvFNPR
mvNNI
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvSubtreeSwap(topology, weight=taxa.size()) )
## references
- citation: Höhna S & Drummond AJ (2012). Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Systematic Biology, 61(1):1-11.
  doi: 10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-abstract/61/1/1/1676649
