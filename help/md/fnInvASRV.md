## name
fnInvASRV
## title
fnInvASRV
## description
Add an invariable-sites component to a submodel.
## details
This model specifies that some fraction pInv of sites are invariable.
If the submodel has m components, this function will return a model with
m+1 components.
## authors
Benjamin Redelings
## see_also
fnUnitMixture
fnGammaASRV
fnMixtureASRV
## example
        p_inv ~ dnUniform(0,1)
        M := fnInv( fnJC(4), p_inv)
        M := fnJC(4) |> fnInv(p_inv)  # Nested functions can be expressed using pipes.

        M := fnJC(4) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
        M := fnJC(4) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (4*2) components - slower.

        # Not recommended -- illustration only.  3 components.
        M := fnJC(4) |> fnInv(p2) |> fnInv(p2) # Fraction of invariable sites is p2 + (1-p2)*p2
## references
