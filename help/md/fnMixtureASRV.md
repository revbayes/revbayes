## name
fnMixtureASRV
## title
fnMixtureASRV
## description
Constructs a mixture model from a collection of submodels.
## details
Each site will evolve according to one of the submodels. The probability that
each site follows a particular submodel is specified by the fractions parameter.
If the rates parameter is given, each submodel rate is multipled by the corresponding
element of the rates vector.

The number of components is the sum of the components of the individual submodels.

## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnGammaASRV
fnInvASRV
## example
        # Two components with different frequencies
        pi1 ~ dnDirichlet([1,1,1,1])
        pi2 ~ dnDirichlet([1,1,1,1])
        weights ~ dnDirichlet([1,1])
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights)
        seq ~ dnPhyloCTMC(psi, M, type="DNA")

        # A free-rates model
        Q := fnJC(4)
        rates ~ dnDirichlet([1,1,1,1])
        weights ~ dnDirichlet([2,2,2,2])
        M := fnMixtureASRV([Q,Q,Q,Q], weights, rates*4)
        
        # Adding rate variation to the frequency-variation model.
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)

## references
