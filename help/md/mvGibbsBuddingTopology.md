## name
mvGibbsBuddingTopology

## title
Gibbs draw of a budding topology

## description
Redraws the whole budding (asymmetric speciation) topology of a `dnFBDSP` tree, holding the species ranges fixed.

## details
Each lineage is assigned an ancestor drawn uniformly from those alive at its birth. Conditional on the ranges every compatible tree carries the same probability, which is what the range process `dnFBDRP` states as a factor of gamma per taxon, so the draw is from the exact conditional. The move is a Gibbs step: the ratio is one and the proposal is always accepted.

The other topology moves rearrange one branch at a time and reach these trees by random walk; this one lands anywhere in the compatible set in a single step. It also covers which lineage continues the ancestral species, since naming an ancestor fixes that, so `mvRotateNode` is redundant beside it.

Ranges and node ages are untouched, so the move needs company: something to sample the ages, such as `mvFossilTipTimeUniform` for the extinction times and `mvNodeTimeSlideUniform` with `mvRootTimeSlideUniform` for the speciation times.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvRotateNode
mvFNPR

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)
rec ~ dnFossilRecord(ranges=tr, complete=false)
rec.clamp(taxa)

moves.append( mvGibbsBuddingTopology(tr, weight=taxa.size()) )
moves.append( mvFossilTipTimeUniform(tr, weight=taxa.size()) )
moves.append( mvNodeTimeSlideUniform(tr, weight=taxa.size()) )
moves.append( mvRootTimeSlideUniform(tr, origin, weight=2) )

## references
