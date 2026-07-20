## name
mvWeightedSubtreeSwap
## title
Weighted Subtree Swap move.
## description
Tree topology move that performs a Weighted Subtree Swap on an unrooted tree.
## details
The basic idea of the `mvWeightedSubtreeSwap` is to propose new trees based on the marginal probabilities (marginal in the sense of integrating over branch lengths for the affected branches). It may be possible that some topology moves are reject only because of the specific branch lengths but otherwise the new topology should have a higher posterior probability.

`mvWeightedSubtreeSwap` randomly picks an internal nodes (clade). Then it finds all possible other clades it can be swap with. Now it tries all these possible swaps and computes the resulting joint probability density. Importantly, it integrates over all possible fractions of the branch length for the swapped node. This is done by breaking the fraction into intervals which are computed by the quantiles of a Beta distribution. The integral is computed by using a piecewise linear extrapolation (similar to a Riemann integral).
## authors
Sebastian Höhna
## see_also
mvFNPR
mvNNI
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    topology ~ dnUniformTopology(taxa)
    moves.append( mvWeightedSubtreeSwap(topology, weight=taxa.size()) )

## references
- citation: Höhna S & Drummond AJ (2012). Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Systematic Biology, 61(1):1-11.
  doi: 10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-abstract/61/1/1/1676649
