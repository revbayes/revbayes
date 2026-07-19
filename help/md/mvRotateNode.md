## name
mvRotateNode

## title
Node rotation move

## description
Permutes the children of a random internal node, leaving every node age and the clade set unchanged.

## details
Child order carries no meaning in an unlabelled tree, so under most distributions this proposes the state the chain is already in and always accepts. The move is for a process that reads the order as state.

`dnFBDSP` is one: it takes a node's first child as the lineage that continues its ancestor's species and the rest as budding (asymmetric speciation) descendants, so rotating a node is a different budding history over the same topology, with different birth times and a different probability. The topology moves reach those states only as a side effect of rearranging the tree, and never propose one on its own.

The permutation is drawn uniformly over the orders that differ from the current one, in both directions, so the proposal is symmetric and carries no Hastings ratio. On a bifurcating node it exchanges the two children. A node with a single child is never chosen.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvNNI
mvFNPR

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)

# which lineage continues the ancestral species at each speciation event
moves.append( mvRotateNode(tr, weight=taxa.size()/2) )

## references
