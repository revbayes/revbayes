## name
fnGammaASRV
## title
fnGammaASRV
## description
Add Gamma-distributed across-site rate variation (ASRV) to a site model.

## details
Each site evolves according to the specified site model, but at an unknown rate
that is Gamma distributed. If the site model parameter is a mixture model with
m components, this function will return a mixture with m*n components.

The continuous Gamma distribution is approximated with a mixture distribution
over n discrete rates, each with probability 1/n.  The Gamma distribution is
constrained to have a mean of 1, so as not to change the  branch lengths.
It therefore has only a single parameter alpha -- the shape parameter.
        - As alpha approaches infinity, all rates across sites become equal (rate variation goes to 0).
        - If alpha = 1, then the rate is exponentially distributed.  Rate variation is substantial.
        - As alpha approaches zero, many sites have rate 0, and many sites have a high rate.

RateMatrix and RateGenerator site model parameters will automatically be converted to a
SiteMixtureModel with a single component.

## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnInvASRV
fnScale
fnMixtureASRV
fnDiscretizeGamma
## example
        # fnGammaASRV( ) constructs a mixture model that represents both the underlying
        #   rate matrix and Gamma-distributed rate variation.
        for (i in 1:10) { taxa[i] = taxon("T"+i) }
        psi ~ dnBDP(lambda=1, rootAge=1, taxa=taxa)
        alpha ~ dnExp(1/10)
        er ~ dnDirichlet( [1,1,1,1,1,1] )
        pi ~ dnDirichlet( [1,1,1,1] )
        M := fnGammaASRV( fnGTR(er, pi), alpha, 4)
        seq ~ dnPhyloCTMC(psi, M, type="DNA",nSites=10)

        # As an alternative approach, models can be built up iteratively using pipes.
        M := fnGTR(er,pi) |> fnGammaASRV(alpha, 4)  

        M := fnGTR(er,pi) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
        M := fnGTR(er,pi) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (2*4) components - slower.

        # The site model parameter can be a mixture model
        weights ~ dnDirichlet([1,1])
        pi1 ~ dnDirichlet( [1,1,1,1,1,1 ] )
        pi2 ~ dnDirichlet( [1,1,1,1,1,1 ] )
        M := fnMixtureASRV([fnGTR(er,pi1),fnGTR(er,pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)

## references
- citation: Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA sequences with variable rates
      over sites: approximate methods
  doi: https://doi.org/10.1007/BF00160154
