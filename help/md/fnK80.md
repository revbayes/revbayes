## name
fnK80

## title
The Kimura (1980) nucleotide rate matrix

## description
DNA evolution model proposed in Kimura (1980).

## details
In this model, all nucleotides have an equal stationary frequency, and transition and transversion rates are allowed to be different. Its only parameter, kappa, codes for the ratio between the rate of transitions and transversions.

The K80 rate matrix elements will be of the form:
    Q[i, j] = c * kappa, if i<->j is a transition
            = c, if i<->j is a transversion

where c is a constant needed to normalize the average rate to 1.

## authors
## see_also
fnJC
fnF81
fnK81
fnT92
fnHKY
fnGTR

## example
    # the ratio between rates of transitions and transversions
    kappa ~ dnExp(0.5)

    # create a K80 rate matrix
    Q := fnK80(kappa)

## references
- citation: Kimura M (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences". Journal of Molecular Evolution. 16:111–20.
  doi: https://doi.org/10.1007/BF01731581
  url: https://link.springer.com/article/10.1007/BF01731581
