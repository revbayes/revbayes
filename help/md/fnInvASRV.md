## name
fnInvASRV
## title
fnInvASRV
## description
Add an invariable-sites component to a site model.
## details
This model specifies that some fraction pInv of sites are invariable.
If the site model parameter is a mixture model with m components, this function will return a model with
m+1 components.
## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnGammaASRV
fnMixtureASRV
## example
        # fnInvASRV( ) creates a mixture model by adding invariant sites to an underlying site model.
        p_inv ~ dnUniform(0,1)
        M := fnInv( fnJC(4), p_inv)
        seq ~ dnPhyloCTMC(psi, M, type="DNA")

        # As an alternative approach, models can be built up iteratively using pipes.
        M := fnJC(4) |> fnInv(p_inv)

        M := fnJC(4) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
        M := fnJC(4) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (4*2) components - slower.

        # Not recommended -- illustration only.  3 components.
        M := fnJC(4) |> fnInv(p2) |> fnInv(p2) # Fraction of invariable sites is p2 + (1-p2)*p2
## references