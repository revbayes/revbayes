## name
SiteMixtureModel
## title
SiteMixtureModel
## description
A weighted collection of discrete character evolution models.
## details
The SiteMixtureModel datatype is a mixture distribution where each
component is a model of discrete character evolution.  Each character evolves
according to one of the component models.  However, the specific model for each
character is not specified in advance.  Instead, each character has some
probability of choosing each component.  These probabilities are specified by
the mixture weights.

## authors
Ben Redelings
## see_also
## example
        M := fnInvASRV(fnGammaASRV(fnJC(4),alpha=1),pInv=0.1)
        M.weights()
        M.nComponents()
        M.rootFrequencies(1)

        # It possible to express nested models using pipes.
        M := fnJC(4) |> fnGammaASRV(alpha=1) |> fnInvASRV(pInv=0.1)
## references
