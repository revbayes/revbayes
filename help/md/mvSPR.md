## name
mvSPR
## title
Subtree Prune and Regraft (SPR) move.
## description
Tree topology move that performs a Subtree Prune and Regraft (SPR) on
an unrooted tree.
## details
`mvSPR` changes tree topology by cutting off a subtree and reattaching it
elsewhere in the original tree using the same subtree branch that was
originally cut. Every unrooted tree of n taxa has 2(n - 3)(2n - 7) SPR
"neighbors" that are one SPR move away (Allen & Steel 2001). This neighborhood
is larger than, and inclusive of, the neighborhood induced by Nearest-Neighbor
Interchange (`mvNNI`). As a result, `mvSPR` is more computationally demanding
than `mvNNI` and may exhibit lower acceptance rates, but explores a broader
range of different topologies and is less likely to get stuck in local optima.
The `mvSPR` move can be applied only to `BranchLengthTreee` objects.
An analogous move for `TimeTree` objects (Fixed Node-height Prune and Regraft;
FNPR) is implemented in `mvFNPR`.
## authors
## see_also
mvFNPR
mvNNI
mvSubtreeSwap
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    moves = VectorMoves()
    
    topology ~ dnUniformTopology(taxa)
    moves.append( mvSPR(topology, weight=taxa.size()) )

## references
- citation: Allen BL, Steel M (2001). Subtree transfer operations and their induced metrics on evolutionary trees. Annals of Combinatorics, 5:1-15.
  doi: 10.1007/s00026-001-8006-8
  url: https://link.springer.com/article/10.1007/s00026-001-8006-8
- citation: Swofford DL, Olsen GJ (1990). Phylogeny reconstruction. Pp. 411–501 in Hillis DM, Moritz C, eds. Molecular Systematics, 1st ed. Sunderland, MA: Sinauer Associates.


