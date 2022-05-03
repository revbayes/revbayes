## name
dnFossilizedBirthDeathRangeMatrix
## title
## description
The fossilized birth death range matrix process (FBDRMatrix) describes the distribution of a matrix of species origination and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. If only a subset of fossil occurrences is provided for each species (e.g. only first and last occurrencces), then the remaining number of fossil samples may be marginalized by specifying `complete=FALSE`. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which are automatically marginalized during when the model is sampled using MCMC. To disable this behavior, use `resample=FALSE`.
## authors
Walker Pett
## see_also
dnFossilizedBirthDeathRange
dnBirthDeathSamplingTreatment
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRMatrix(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, complete=FALSE)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
