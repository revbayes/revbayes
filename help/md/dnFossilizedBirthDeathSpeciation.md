## name
dnFossilizedBirthDeathSpeciation
## title
## description
The fossilized birth death speciation process (FBDSP) describes the diversification and sampling of extant and extinct species trees under a mixed model of asymmetric, symmetric and anagenetic speciation.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. If only a subset of fossil occurrences is provided for each species (e.g. only first and last occurrencces), then the remaining number of fossil samples may be marginalized by specifying `complete=FALSE`. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleFBD`. Tips represent extinction events, and therefore should be sampled during MCMC using e.g. `mvTipTimeSlideUniform`.
## authors
Walker Pett
## see_also
dnFossilizedBirthDeathRange
dnBirthDeathSamplingTreatment
mvResampleFBD
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

origin ~ dnExp(1/10)

bd ~ dnFBDSP(originAge=min_age+origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, complete=FALSE)

moves.append( mvFNPR(bd, weight = taxa.size()) )
moves.append( mvNodeTimeSlideUniform(bd, weight = taxa.size()) )
moves.append( mvRootTimeSlideUniform(bd, origin=origin, weight = taxa.size()) )
moves.append( mvTipTimeSlideUniform(bd, weight = taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
