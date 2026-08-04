## name
dnBDS
## title
## description
The birth-death-with-rateshifts (BDS) process of Silvestro et al. (2019): a fossilized birth-death range process under complete lineage sampling.
## details
This is the birth-death range process on its own, like `dnFBDRP`, but under the assumption of complete lineage sampling: lineages are treated as independent, so there is no coexistence (gamma) factor, and each range is normalized by the fossil non-detection probability over its interval. It describes how species ranges diversify but does not include the probability of the fossil record itself; pair it with a `dnFossilRecord` node clamped to the occurrences to get the full model, choosing `complete=TRUE` or `complete=FALSE` (first/last) there.

Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleAugmentedAges`. By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

The constructor is also available under the alias `dnBirthDeathWithRateshifts`.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
mvStratigraphicRange
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

# the birth-death-with-rateshifts range process
bd ~ dnBDS(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, complete=false)
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.
	  doi: https://doi.org/10.1017/pab.2019.23
	  url: https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C
