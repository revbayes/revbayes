## name
dnFossilizedBirthDeathSpeciation
## title
The fossilized birth-death range process over trees
## description
The fossilized birth death speciation process (FBDSP) describes the diversification and sampling of extant and extinct species trees under a mixed model of asymmetric (budding) and anagenetic speciation.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Like `dnFBDRP`, this is the birth-death range process on its own (here over trees); pair it with a `dnFossilRecord` node clamped to the occurrences to add the probability of the fossil record, choosing `complete=TRUE` or `complete=FALSE` (first/last) there. Under the hood, the fossil data is augmented with the first and last appearance of each species. A tree carries divergence times only, so these are state of the distribution rather than elements of its value and no tree move reaches them: sample them with `mvStratigraphicRange`. By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

Every age in the tree needs a move, and one left without a move is silently held at its initial value. Tips are extinction events and may fall below their taxon's occurrence range; `mvFossilTipTimeUniform` reads that off the tree and samples them between the present and the youngest occurrence. The root age is the first speciation event and is sampled like any other node age, but `mvNodeTimeSlideUniform` skips the root, so pair it with `mvRootTimeSlideUniform`.
Symmetric (bifurcating) speciation is not available: the configurations it allows have no move, so a positive symmetric speciation probability would leave part of the state space unreachable.

## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
mvFossilTipTimeUniform
mvRootTimeSlideUniform
mvStratigraphicRange
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

min_age = 0.0
for(i in 1:taxa.size())
{
	if ( taxa[i].getMinAge() > min_age )
	{
		min_age = taxa[i].getMinAge()
	}
}

origin_offset ~ dnExp(1/10)
origin := min_age + origin_offset
moves.append( mvSlide(origin_offset, weight = 2) )

bd ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, conditioned on the range process and clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, complete=false)
rec.clamp(taxa)

moves.append( mvFNPR(bd, weight = taxa.size()) )
moves.append( mvNodeTimeSlideUniform(bd, weight = taxa.size()) )
moves.append( mvRootTimeSlideUniform(bd, origin=origin, weight = taxa.size()) )
moves.append( mvFossilTipTimeUniform(bd, weight = taxa.size()) )

# which lineage continues its ancestor's species at each speciation event
moves.append( mvRotateNode(bd, weight = taxa.size()) )

# the first and last appearances, which the tree does not carry
moves.append( mvStratigraphicRange(bd, weight = taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
