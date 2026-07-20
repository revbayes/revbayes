## name
mvGibbsSubtreeSwap
## title
Gibbs Subtree Swap move.
## description
Tree topology move that performs a Gibbs Subtree Swap on an unrooted tree.
## details
`mvGibbsSubtreeSwap` randomly picks a first internal node (clade). Then it finds all possible other clades it can be swap with. Now it tries all these possible swaps and computes the resulting joint probability density. The final proposal is done proportional to the probability density.
The move is thus similar to the basic subtree swap, which randomly pick the second clade.
This move can be efficient to reach quicker the stationary distribution, but is also more costly in computational time.
## authors
Sebastian Höhna
## see_also
mvFNPR
mvNNI
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvGibbsSubtreeSwap(topology, weight=taxa.size()) )
## references
- citation: Höhna S & Drummond AJ (2012). Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Systematic Biology, 61(1):1-11.
  doi: 10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-abstract/61/1/1/1676649
