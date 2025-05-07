## name
fnDiscretizeDistribution
## title
Discretize a Continuous Distribution
## description
`fnDiscretizeDistribution` transforms a continuous probability distribution into a discrete one by dividing it into a specified number of categories.

## details
This function takes as two arguments: a continuous probability distribution and a specified number of categories (`num_cats`).
It then yields a sequence of median values that approximate the distribution, assuming that each bin represents an equal probability mass.

## authors
## see_also
fnDiscretizeGamma

fnDiscretizeBeta

## example
    # Using a Normal distribution to discretize it into 5 categories
    discrete_values := fnDiscretizeDistribution( dnNormal( 0.0, 1.0 ), 5 )

    # print the discretized values to the screen
    discrete_values

## references
