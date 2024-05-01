## name
fnUnitMixture
## title
fnUnitMixture
## description
Create a SiteMixtureModel from a RateMatrix or RateGenerator
## details
This function creates a SiteMixtureModel with one component by specifying the
rate and root frequencies for a RateGenerator.  The rate defaults to 1, leaving
the underlying model unchanged.

If the submodel is a RateMatrix, the root frequencies default to the equilibrium
frequencies of the RateMatrix.  However, a RateGenerator may not have equilibrium
frequencies, so the root frequencies must be specified explicitly.

In many cases it is not necessary to explicitly call fnUnitMixture(), RevBayes can
automatically convert a RateMatrix to a SiteMixtureModel.

## authors
## see_also
## example
        M := fnUnitMixture( fnJC(4) )
        M := fnJC(4) |> fnUnitMixture()  # nested functions can be expressed using pipes.
        
        # Explicit conversion to SiteMixtureModel
        M := fnGTR(er,pi) |> fnUnitMixture() |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)
        # Implicit conversion to SiteMixtureModel
        M := fnGTR(er,pi) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)
        
        # Starting the model at non-equilibrium frequencies.
        M := fnDECRateMatrix(dr,er,"Include") |> fnUnitMixture(rootFrequencies=simplex(rep1,n_states))
        
## references
