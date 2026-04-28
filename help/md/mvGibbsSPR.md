## name
mvGibbsSubtreePruneAndRegraft
## title
Gibbs Subtree Prune and Regraft (GPR) move.
## description
Tree topology move that performs a Gibbs Subtree Prune and Regraft (GPR) on
an unrooted tree.
## details
`mvGibbsSubtreePruneAndRegraft` basically performs a SPR but tries all possible reattachment points, evaluates these based on the current joint DAG probability density, and then proposes the reattachment based on the probability density.
This move can be efficient to reach quicker the stationary distribution, but is also more costly in computational time.
## authors
Sebastian Höhna
## see_also
mvGibbsFNPR
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()

    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))
    moves.append( mvGibbsSubtreePruneAndRegraft(topology, weight=taxa.size()) )
## references
- citation: Höhna S & Drummond AJ (2012). Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Systematic Biology, 61(1):1-11.
  doi: 10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-abstract/61/1/1/1676649
