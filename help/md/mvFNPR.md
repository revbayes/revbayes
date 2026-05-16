## name
mvFNPR
## title
Fixed Node-height Prune and Regraft (FNPR) move.
## description
Tree topology move that prunes and re-attaches a subtree without changing any
node heights.
## details
`mvFNPR` randomly picks node i which is neither a tip nor the root, and prunes
the subtree originating with this node. It then picks another node j such that
j is younger than i but the parent of j is older than i, and re-attaches the
pruned subtree onto the branch above j. Because the node height of i is fixed
rather than re-adjusted, the FNPR move represents a special case of the fully
general time tree version of the subtree prune and regraft (SPR) move. This
fully general version is also known as the Wilson-Balding move. `mvFNPR` often
exhibits higher acceptance rates than the Wilson-Balding move.
## authors
Sebastian Höhna
## see_also
mvSPR
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    height ~ dnUniform(0, 10)
    moves = VectorMoves()

    # Simulate a simple TimeTree
    tree ~ dnBDP(lambda=1.0, mu=0.2, rootAge=height, taxa=taxa)

    # Assign it a mvFNPR move
    moves.append( mvFNPR(tree, weight=taxa.size()) )

## references
- citation: Höhna S, Defoin-Platel M, Drummond AJ (2008). Clock-constrained tree proposal operators in Bayesian phylogenetic inference. 1--7 in 8th IEEE International Conference on BioInformatics and BioEngineering (BIBE 2008). Athens, Greece, October 2008.
  doi: 10.1109/BIBE.2008.4696663
  url: https://alexeidrummond.org/assets/publications/2008-hoehna-clock-bibe.pdf
