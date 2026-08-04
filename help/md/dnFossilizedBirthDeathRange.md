## name
dnFossilizedBirthDeathRange
## title
## description
The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species origination and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
This distribution is the birth-death range process on its own: it describes how species ranges diversify, and the fossil occurrences bound those ranges, but it does not include the probability of the fossil record itself. Pair it with a `dnFossilRecord` node clamped to the occurrences to get the full model. That node also selects the reporting model, i.e. how sampled occurrences reach the record: `complete` (all reported) or first/last (only the oldest and youngest). A range process with no `dnFossilRecord` attached is a valid model, but it does not condition on the fossil record at all.

Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleAugmentedAges`. The related `dnBDS` distribution assumes complete lineage sampling instead (the Birth-Death with Rateshifts model of Silvestro et al. 2019). By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

## authors
June Walker
## see_also
dnFossilRecord
dnBDS
mvStratigraphicRange
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, conditioned on the range process and clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, complete=false)
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
	- citation: Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.
	  doi: https://doi.org/10.1017/pab.2019.23
	  url: https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C
