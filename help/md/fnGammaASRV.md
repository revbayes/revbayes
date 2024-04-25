## name
fnGammaASRV
## title
fnGammaASRV
## description
Add Gamma-distributed across-site rate variation (ASRV) to a submodel.

## details
Each site evolves according to the specified submodel, but at an unknown rate
that is Gamma distributed. If the submodel has m components, this function will
return a model with m*n components.

The continuous Gamma distribution is approximated with a mixture distribution
over n discrete rates, each with probability 1/n.  The Gamma distribution is
constrained to have a mean of 1, so as not to change the  branch lengths.
It therefore has only a single parameter alpha -- the shape parameter.
        - As alpha approaches infinity, rate variation goes to 0.
        - If alpha = 1, then the rate is exponentially distributed.  Rate variation is substantial.
        - As alpha approaches zero, many sites have rate 0, and many sites have a high rate.

RateMatrix and RateGenerator submodels will automatically be converted to a
SiteMixtureModel with a single component.

## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnInvASRV
fnMixtureASRV
fnDiscretizeGamma
## example
        alpha ~ dnExp(1/10) # Don't put too much prior belief on high rate variation.
        M := fnGammaASRV( fnGTR(er, pi), alpha, 4)
        M := fnGTR(er,pi) |> fnGammaASRV(alpha, 4)  # Nested functions can be expressed using pipes.
        seq ~ dnPhyloCTMC(psi, M, type="DNA")

        M := fnJC(4) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
        M := fnJC(4) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (4*2) components - slower.

        # The submodel can be a mixture model
        weights ~ dnDirichlet([1,1])
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)

        # Rate is a product of two Gamma-distributed variables.  For illustration only -- NOT recommended.
        M := fnJC(4) |> fnGammaASRV(alpha1,4) |> fnGammaASRV(alpha2,4)   # This has 16 components.

## references
- citation: Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA sequences with variable rates
      over sites: approximate methods
  doi: https://doi.org/10.1007/BF00160154
