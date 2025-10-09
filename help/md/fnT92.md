## name
fnT92

## title
The Tamura (1992) nucleotide rate matrix

## description
DNA evolution model proposed in Tamura (1992).

## details
In this model, A and T have an equal stationary frequency, with G and C frequencies distinct, and transition and transversion rates are allowed to be different. Its first parameter, kappa, codes for the ratio between the rate of transitions and transversions. Its second parameter, gc, codes for the compound frequency of G and C nucleotides.

The T92 rate matrix elements will be of the form:
    Q[i, j] = c * kappa * gc / 2, if i<->j is a transition and j is C or G
            = c * gc / 2, if i<->j is a transversion and j is C or G
            = c * kappa * (1 - gc) / 2, if i<->j is a transition and j is A or T
            = c * (1 - gc) / 2, if i<->j is a transversion and j is A or T

where c is a constant needed to normalize the average rate to 1.

## authors
## see_also
fnJC
fnF81
fnK80
fnK81
fnHKY
fnTrN
fnGTR

## example
    # the ratio between rates of transitions and transversions
    kappa ~ dnExp(0.5)

    # the frequency of G and C nucleotides
    gc ~ dnUnif(0, 1)

    # create a T92 rate matrix
    Q := fnT92(kappa, gc)

## references
- citation: Tamura K (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases". Molecular Biology and Evolution. 9:678–87.
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040752
  url: https://academic.oup.com/mbe/article/9/4/678/1254082
