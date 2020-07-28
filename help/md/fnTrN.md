## name
fnTrN
## title
The Tamura-Nei (1993) nucleotide rate matrix
## description
The Tamura-Nei nucleotide rate matrix.
## details
## authors
## see_also
fnHKY
## example
        kappaAG ~ dnLognormal(0,1)    # The purine transition rate
        kappaCT ~ dnLognormal(0,1)    # The pyrimindine transition rate
        pi ~ dnDirichlet( v(1,1,1,1) )
        Q := fnTrN(kappaAT, kappaCT, ,pi)
## references
- citation: Tamura, K. and M. Nei. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular biology and evolution (1993) 10(3):512-526.
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040023
  url: https://academic.oup.com/mbe/article/10/3/512/1016366
