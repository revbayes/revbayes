## name
mvBuddingTopologyRedraw

## title
Gibbs draw of a budding topology

## description
Redraws the whole budding (asymmetric speciation) topology of a `dnFBDSP` tree, holding the species ranges fixed.

## details
Each lineage is assigned an ancestor drawn uniformly from those alive at its birth. Conditional on the ranges every compatible tree carries the same probability, which is what the range process `dnFBDRP` states as a factor of gamma per taxon, so the draw is from the exact conditional and every proposal is accepted.

That equality needs pure budding. With `lambda_a > 0` an anagenetic attachment sits exactly at the ancestor's extinction time and carries `lambda_a` in place of `lambda` and `mu`, so compatible trees no longer share one density. The uniform draw cannot reach such an attachment, and rebuilding a budding tree would discard every sampled ancestor, so the move refuses to be set up on a process with a positive anagenetic rate. Use the MH topology moves there.

That holds while the tree carries no character data. A phylogenetic likelihood breaks the tie between compatible trees, and the move becomes an independence proposal drawn from the tree prior rather than a Gibbs step. It stays valid, since the draw is uniform over the compatible set and independent of the current tree, so the ratio is still one and the likelihood alone decides acceptance. It stops being efficient: a blind draw rarely agrees with an informative alignment, so pair it with `mvFNPR` or `mvNNI` for the local rearrangements a likelihood rewards.

The other topology moves rearrange one branch at a time and reach these trees by random walk; this one lands anywhere in the set compatible with the current birth and death ages in a single step.

It rebuilds each budding node at the birth age the taxon already has, so it redraws which lineage a species buds from but never which species buds. The oldest birth in particular is taken as the root lineage and is never reattached. Pair it with `mvRotateNode`, which swaps the roles at a node and is the only move that relabels which lineage carries a given birth age. Without it the oldest birth stays on one taxon for the whole run and every per-taxon birth age is sampled too narrowly.

Ranges and node ages are untouched, so the move needs company for those too: `mvFossilTipTimeUniform` for the extinction times and `mvNodeTimeSlideUniform` with `mvRootTimeSlideUniform` for the speciation times.

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

moves.append( mvBuddingTopologyRedraw(tr, weight=taxa.size()) )
# relabels which lineage carries each birth age; the Gibbs move cannot
moves.append( mvRotateNode(tr, weight=taxa.size()) )
moves.append( mvFossilTipTimeUniform(tr, weight=taxa.size()) )
moves.append( mvNodeTimeSlideUniform(tr, weight=taxa.size()) )
moves.append( mvRootTimeSlideUniform(tr, origin, weight=2) )

## references
