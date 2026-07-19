## name
dnFossilizedBirthDeathSpeciation
## title
## description
The fossilized birth death speciation process (FBDSP) describes the diversification and sampling of extant and extinct species trees under a mixed model of asymmetric, symmetric and anagenetic speciation.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Like `dnFBDRP`, this is the birth-death range process on its own (here over trees); pair it with a `dnFossilRecord` node clamped to the occurrences to add the probability of the fossil record, choosing `complete=TRUE` or `complete=FALSE` (first/last) there. Under the hood, the fossil data is augmented with oldest occurrence ages for each species. These ride with the moves on this node while `resample=TRUE`, and may also be sampled explicitly with `mvResampleAugmentedAges`.

Every age in the tree needs a move, and one left without a move is silently held at its initial value. Tips are extinction events and may fall below their taxon's occurrence range, so sample them with `mvExtendedTipTimeUniform` rather than `mvFossilTipTimeUniform`. The root age is the first speciation event and is sampled like any other node age, but `mvNodeTimeSlideUniform` skips the root, so pair it with `mvRootTimeSlideUniform`.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
mvExtendedTipTimeUniform
mvRootTimeSlideUniform
mvResampleAugmentedAges
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
moves.append( mvExtendedTipTimeUniform(bd, weight = taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
