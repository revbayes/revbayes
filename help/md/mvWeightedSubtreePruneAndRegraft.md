## name
mvWeightedSubtreePruneAndRegraft
## title
Weighted Subtree Prune and Regraft move.
## description
Tree topology move that performs a Weighted Subtree Prune and Regraft on an unrooted tree.
## details
The basic idea of the `mvWeightedSubtreePruneAndRegraft` is to propose new trees based on the marginal probabilities (marginal in the sense of integrating over branch lengths for the affected branches) of all possible other trees after the subtree for pruning has been chosen (see `mvGibbsSubtreePruneAndRegraft`). It may be possible that some topology moves are reject only because of the specific branch lengths but otherwise the new topology should have a higher posterior probability.

`mvWeightedSubtreePruneAndRegraft` randomly picks a node and the corresponding branch leading to it (subtree). Then it finds all possible other branches (or nodes, depending on your viewpoint) it can be reattached to. Now it tries all these possible reattachments and computes the resulting joint probability density (see `mvWeightedBranchLengthScale`). Importantly, it integrates over all possible fractions of the branch length it reattaches to. This is done by breaking the fraction into intervals which are computed by the quantiles of a Beta distribution. The integral is computed by using a piecewise linear extrapolation (similar to a Riemann integral).
## authors
Sebastian Höhna
## see_also
mvWeightedBranchLengthScale
mvGibbsSubtreePruneAndRegraft
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvWeightedSubtreePruneAndRegraft(topology, weight=taxa.size()) )
## references
