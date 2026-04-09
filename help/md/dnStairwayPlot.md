## name
dnStairwayPlot
## title
StairwayPlot Distribution
## description
Bayesian StairwayPlot for inferring single population demographic histories
from site frequency spectra.
## details
The `StairwayPlot Distribution` specifies a distribution on the site frequency
spectrum of a single panmictic population. The site frequency spectrum varies
depending on the demographic history, that is, changes in the effective
population size. You can pass in a vector of effective population sizes, which
are assumed to change exactly at the expected coalescent events (see the
references for the description of the theory). The name of the distribution
comes from the stairway-like stepwise function of the effective population
sizes.

The most important arguments are:
`theta`: The vector of effective population sizes in units of 4 * Ne * mu. This
    vector must be of size N-1, where N is the number of individuals.
`numSites`: The number of sites used. This should be the same as summing the
    SFS, but is required for simulation/initialization.
`numIndividuals`: The number of individuals, which corresponds to the number
    of bins in the SFS. The fixed sites are not considered. This should be the
    same as in the observed SFS, but is required for simulation/initialization.
`folded`: Whether the likelihood is computed based on a folded SFS or not.
`monomorphicProbability`: How we should compute the probability of monomorphic
    sites. See the references for details.
`coding`: Whether we are conditioning on not including monomorphic sites
    ("no-monomorphic"), not including singletons ("no-singletons"), or whether
    we have "all" sites.

## authors
Sebastian Höhna
## see_also
## example
	# let's assume we have some SFS "observed"
	obs_sfs = [ 305082, 44248, 32223, 28733, 28220, 26205, 27477, 26618, 27533,
                26945, 28736, 28671, 31277, 31250, 34352, 34859, 38331, 40005,
                45666, 48986, 65829, 64363, 70895, 74114, 82705, 88226, 102194,
                114566, 130176, 143775, 169216, 191624, 230016, 276489, 333069,
                394810, 501961, 653809, 890077, 1349350, 50296796 ]

	# we need to remove the fixed sites
	obs_sfs[obs_sfs.size()] <- 0

	# get the number of individuals and the number of sites
	N_IND   = abs(obs_sfs.size()-1)
	N_SITES = round(sum(obs_sfs))

	# now specify a different theta per interval
	for (i in 1:(N_IND-1)) {
	  theta[i] ~ dnUnif(0.0, 0.1)
	}

	sfs ~ dnStairwayPlot( theta, numSites=N_SITES, numIndividuals=N_IND, folded=TRUE,
                          monomorphicProbability="rest", coding="all" )
	sfs.clamp( obs_sfs )

## references
- citation: Liu X, Fu Y-X (2015). Exploring population size changes using SNP frequency spectra. Nature Genetics, 47:555--559.
  doi: 10.1038/ng.3254
  url: https://pmc.ncbi.nlm.nih.gov/articles/PMC4414822/pdf/nihms-668186.pdf
- citation: Höhna S, Catalán A (2025). Bayesian StairwayPlot for inferring single population demographic histories from site frequency spectra. Molecular Ecology Resources, 25:e14087.
  doi: 10.1111/1755-0998.14087
  url: https://onlinelibrary.wiley.com/doi/pdf/10.1111/1755-0998.14087
