## name
fnTrN

## title
The Tamura-Nei (1993) nucleotide rate matrix

## description
DNA evolution model proposed in Tamura & Nei (1993).

## details
In this model, nucleotide base frequencies are different, and the two transition rates (A <-> G and C<->T) can be different to each other, and to the transversion rate. The first argument, kappa1, defines the ratio between the rate of A <-> G (i.e. purine) transitions to transversions. The second argument, kappa2, defines the ratio between the rate of C <-> T (i.e. pyrimidine) transitions to transversions. The third argument, baseFrequencies, defines the stationary frequencies of nucleotide bases. 

The TrN rate matrix elements are of the form:
    Q[i, j] = c * kappa1 * baseFrequencies[j], if i<->j is A<->G
            = c * kappa2 * baseFrequencies[j], if i<->j is C<->T
            = c * baseFrequencies[j], otherwise

where c is a constant needed to normalize the average rate to 1 

## authors
## see_also
fnJC
fnK80
fnK81
fnT92
fnHKY
fnGTR

## example
    # A <-> G transition rate
    kappaAG ~ dnLognormal(0,1)

    # C <-> T transition rate
    kappaCT ~ dnLognormal(0,1)

    # nucleotide base frequencies
    pi ~ dnDirichlet( v(1,1,1,1) )

    # create a TrN rate matrix
    Q := fnTrN(kappaAG, kappaCT, ,pi)

## references
- citation: Tamura, K. and M. Nei (1993). "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees". Molecular biology and evolution. 10(3):512-526.
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040023
  url: https://academic.oup.com/mbe/article/10/3/512/1016366
