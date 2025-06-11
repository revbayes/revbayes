## name
fnScale
## title
fnScale
## description
Scale a vector of SiteMixtureModels
## details
This function has two forms.  The first form takes a SiteMixtureModel `model` and scales it by
a rate `rate`.  This form returns SiteMixtureModel.

The second form takes SiteMixtureModel[] `models` and RealPos[] `rates`, and scales `models[i]`
by `rates[i]`.  This form returns SiteMixtureModel[].

As a shortcut, if the second argument `rates` is a vector but the first element `model` is not,
then the first argument will be automatically replaced with a vector of SiteMixtureModels of the
same length as `rates`, where each element is identical to `model`.

## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnInvASRV
fnMixtureASRV

## example
        Q = fnJC(4)                    # The rate of Q is 1

        # Operating on SiteMixtureModel
        Q2 = fnScale(Q,2)              # The rate of Q2 is 2

        # Operating on SiteMixtureModel[]
        Qs = fnScale([Q,Q],[1,2])      # Qs[1] and Qs[2] have rates 1 and 2
        Qs = fnScale(Q,    [1,2])      # An abbreviation for the above.

        # We can build up models iteratively using pipes
        Qs = Q |> fnScale([1,2])       # A shorter abbreviation.

        # A JC+LogNormal[4] ASRV model
        site_rates := dnLognormal(0,lsigma) |> fnDiscretizeDistribution(4)
        MM := fnJC(4) |> fnScale(site_rates) |> fnMixtureASRV()
        M := fnScale(MM, 1/MM.rate())

        # A FreeRates[5] ASRV model
        rates ~ dnDirichlet( [1,1,1,1,1] )
        weights ~ dnDirichlet( [2,2,2,2,2] )
        MM := fnJC(4) |> fnScale(rates) |> fnMixtureASRV(weights)
        M := fnScale(MM, 1/MM.rate())
## references
