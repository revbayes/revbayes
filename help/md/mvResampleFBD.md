## name
mvResampleFBD
## title
## description
This move resamples an oldest occurrence age for a random species in a fossilized birth death process described by `dnFBDRP` or `dnFBDRMatrix`
## details
Under the hood, FBD fossil data is augmented with oldest occurrence ages for each species, which are automatically marginalized during when the model is sampled using MCMC. These ages can also be resampled manually using this move.
## authors
Walker Pett
## see_also
dnFossilizedBirthDeathRange
dnFossilizedBirthDeathRangeMatrix
## example
bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, resample=FALSE)

moves.append( mvResampleFBD(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
