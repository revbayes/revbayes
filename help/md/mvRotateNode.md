## name
mvRotateNode

## title
Species continuation move

## description
Moves the continuation of an ancestral species between the children of a random internal node, leaving every node age and the topology unchanged.

## details
Under budding (asymmetric) speciation exactly one child of each node continues its ancestor's species, and the others begin new species there. Which one continues is a free parameter, separate from the topology and from the node ages, and this is the move that samples it. It is recorded on the child rather than in the child order, so no topology move can reassign it silently.

A node with no continuing child, or with more than one, is not a state `dnFBDSP` can produce, so the move acts only on nodes that currently name exactly one, and always leaves exactly one. The replacement is drawn uniformly among the remaining children, and the reverse draw is from a set of the same size, so the proposal is symmetric and carries no Hastings ratio.

The continuation determines each range's origination time, so the topology moves reach these states only as a side effect of rearranging the tree, and never propose one on its own.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvBuddingTopologyRedraw
mvStratigraphicRange

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)

# which lineage continues the ancestral species at each speciation event
moves.append( mvRotateNode(tr, weight=taxa.size()/2) )

## references
