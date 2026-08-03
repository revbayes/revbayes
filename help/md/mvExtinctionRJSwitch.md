## name
mvExtinctionRJSwitch
## title
Jump a taxon's extinction time to and from the present
## description
Moves one taxon's extinction time between the present, where it survived to the present without
being sampled, and a time above it.
## details
When rho is below one, a taxon reported extinct may have survived unsampled. Its extinction time
then has a point mass at the present, carrying weight 1 - rho, mixed with a density over times
above it. The element moves are continuous, so the point has probability zero under them and is
never proposed. This move jumps between the two states, which is what makes rho < 1 samplable.

One direction is deterministic and the other draws uniformly below the taxon's youngest
appearance, so the Hastings ratio is the width of that support. Taxa reported extant are pinned at
the present and are never chosen.

At rho = 1 the point mass has no weight and the move is wasted effort, not an error.
## authors
## see_also
mvMatrixElementSlide
mvStratigraphicRange
mvRJSwitch
## example
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=0.5, timeline=timeline)

    moves.append( mvMatrixElementSlide(bd, delta=1, weight=taxa.size()/10) )
    moves.append( mvExtinctionRJSwitch(bd, weight=taxa.size()/10) )
## references
