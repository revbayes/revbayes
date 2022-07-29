## name
dnFossilizedBirthDeathRange
## title
## description
The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species origination and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. If only a subset of fossil occurrences is provided for each species (e.g. only first and last occurrencces), then the remaining number of fossil samples may be marginalized by specifying `complete=FALSE`. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleFBD`. Setting `BDS` to true causes the model to assume complete lineage sampling i.e. using the Birth-Death with Rateshifts (BDS) model of Silvestro et al. (2019).
## authors
Walker Pett
## see_also
dnBirthDeathSamplingTreatment
mvResampleFBD
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
	- citation: Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.
	  doi: https://doi.org/10.1017/pab.2019.23
	  url: https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C
