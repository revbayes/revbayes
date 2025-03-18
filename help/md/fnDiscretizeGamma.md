## name
fnDiscretizeGamma
## title
Discretize a Gamma distribution
## description
`fnDiscretizeGamma` approximates a continuous gamma distribution by dividing it into a specified number of discrete categories (quantiles), using either the mean or the median of each interval.

## details
This function takes a gamma distribution parameterized by `shape` and `rate`, along with a specified number of categories (`numCats`).
It then discretizes the distribution into `numCats` bins, ensuring that each bin represents an equal probability mass.
The representative values for each category can be chosen based on either the mean or the median of the interval.

## authors
## see_also
fnDiscretizeDistribution

fnDiscretizeBeta

## example
    # to obtain the mean of the quantiles
    alpha = 0.5

    discrete_values := fnDiscretizeGamma( shape = alpha, rate = alpha, numCats = 4, median = FALSE )


## references
- citation: Yang, Z (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: Approximate methods. J Mol Evol 39, 306–314.
  doi: 10.1007/BF00160154
  usl: https://doi.org/10.1007/BF00160154
