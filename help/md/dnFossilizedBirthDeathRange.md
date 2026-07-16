## name
dnFossilizedBirthDeathRange
## title
## description
The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species origination and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
This distribution is the birth-death range skeleton on its own: it describes how species ranges diversify, and the fossil occurrences bound those ranges, but it does not include the probability of the fossil record itself. Pair it with a `dnFossilRecord` node clamped to the occurrences to get the full model. That node also selects the reporting model, i.e. how sampled occurrences reach the record: `complete` (all reported), `firstlast` (only the oldest and youngest), or `uniform` (an exchangeable subset). A skeleton with no `dnFossilRecord` attached is a valid model, but it does not condition on the fossil record at all.

Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleFBDR`. Setting `BDS` to true causes the model to assume complete lineage sampling i.e. using the Birth-Death with Rateshifts (BDS) model of Silvestro et al. (2019).

The deprecated `dnFBDRMatrix` is the older fused form of this process: it carries the fossil-record term internally and takes the occurrences as a constructor argument rather than as clamped data. Its `complete` and `reporting` arguments are replaced by `reporting` on `dnFossilRecord`, where `complete=TRUE` becomes `reporting="complete"`.
## authors
Walker Pett
## see_also
dnFossilRecord
dnBirthDeathSamplingTreatment
mvResampleFBDR
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, conditioned on the skeleton and clamped to the observed occurrences
rec ~ dnFossilRecord(skeleton=bd, reporting="firstlast", taxa=taxa)
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
