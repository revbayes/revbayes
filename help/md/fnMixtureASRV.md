## name
fnMixtureASRV
## title
fnMixtureASRV
## description
Constructs a mixture model from a collection of site models.
## details
Each site will evolve according to one of the input site models, which may also
be mixture models.  The probability that each site follows a particular site model
is specified by the fractions parameter.

The number of components in the resulting mixture model is the sum of the number
of components of the input mixture models.

If the fractions parameter is missing, then each of the given models is given equal
weight.

## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnGammaASRV
fnInvASRV
fnScale

## example
        # Two components with different frequencies
        for (i in 1:10) { taxa[i] = taxon("T"+i) }
        psi ~ dnBDP(lambda=1, rootAge=1, taxa=taxa)
        pi1 ~ dnDirichlet([1,1,1,1])
        pi2 ~ dnDirichlet([1,1,1,1])
        weights ~ dnDirichlet([1,1])
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights)
        seq ~ dnPhyloCTMC(psi, M, type="DNA", nSites=10)

        # A weight of 1/2 on each model because the weights are missing.
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)])

        # Adding rate variation to the frequency-variation model.
        M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)


        
## references
