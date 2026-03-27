## name
mvWeightedBranchLengthScale
## title
Weighted Branch Length Scale move.
## description
Branch Length Scale move that proposes a new branch reattachment based on numeric integration on an unrooted tree.
## details
The basic idea of the `mvWeightedBranchLengthScale` is to move an internal node left or right (or up and down, if the tree is considered vertically). Consider a random node. This node has a branch that leads to it (yes, for unrooted trees we still internally give it a direction for easy identification). Now this node is attached to the remaining tree by that branch (connected through an internal node). So that internal node has two more branches connected to it, we arbitrarily call them left and right (or up and down, if you wish). Imagine you could slide this node left or right, which would change the fraction of the left branch to the right branch, but keeping the total the same. The `mvWeightedBranchLengthScale` move tries to propose new values from an approximated probability density of all left/right fractions. This is done by breaking the fraction into intervals which are computed by the quantiles of a Beta distribution. The integral is computed by using a piecewise linear extrapolation (similar to a Riemann integral).
## authors
Sebastian Höhna
## see_also
mvBranchLengthScale
mvNNI
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvWeightedBranchLengthScale(topology, weight=taxa.size()) )
## references
