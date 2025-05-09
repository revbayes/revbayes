## name
fnF81

## title
The Felsenstein (1981) rate matrix

## description
DNA evolution model proposed in Felsenstein (1981).

## details
In this model, states are allowed to have different stationary frequencies, and exchangeability rates between states are equal. Its only argument, baseFrequencies, codes for said stationary frequencies. While this is usually used for DNA (and therefore has four states), the function can take any number of states, and therefore be used for many other applications (such as aminoacid or morphological evolution).

The F81 rate matrix elements will be of the form:
    Q[i, j] = c * baseFrequencies[j]

where c is a constant needed to normalize the average rate to 1

## authors
## see_also
fnJC
fnK80
fnK81
fnT92
fnHKY
fnTrN
fnGTR

## example
    # stationary base frequencies
    baseFrequencies ~ dnDirichlet(v(1,1,1,1))

    # create an F81 rate matrix
    Q := fnF81(baseFrequencies)

## references
- citation: Felsenstein J (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach". Journal of Molecular Evolution. 17:368–76.
  doi: https://doi.org/10.1007/BF01734359
  url: https://link.springer.com/article/10.1007/BF01734359
