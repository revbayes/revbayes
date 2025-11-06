## name
fnDiscretizeBeta
## title
Disctetize a beta distribution
## description
Select representative values from `num_cats` discrete subdivisions of a beta distribution.

## details

A beta distribution is defined by two shape parameters, alpha and beta.

Where a parameter or prior is defined based on the beta distribution, it may be more tractable to evaluate likelihoods at a fixed number of points from the distribution.  These representative points can be computed using `dnDiscretizeBeta`.

In practice, these values are computed as follows:

Let _n_ be the number of categories.
If `median = TRUE`, the quantile function is performed at the midpoint of each category.  Call this vector _q_.
_q_ is then normalized by dividing against its sum, so its elements sum to one; then multiplied by a factor _n_ * _alpha) / (_alpha_ + _beta_).

The computation to obtain the mean for each category, when `median = FALSE`, is more complex, making use of the incomplete beta function ( Majumder & Bhattacharjee 1973).

A real-world use case is available in Wright et al. (2016), with discussion of the properties of the beta distribution. Corresponding tutorials are available at https://www.palass.org/sites/default/files/media/publications/newsletters/number_106/number_106_0.pdf and https://revbayes.github.io/tutorials/morph_tree/V2.html.


## authors
## see_also
A translation of `fnDiscretizeBeta` into R is available at https://gist.github.com/ms609/883632d10d4d80ea5391cee9c47071fc.

## example
    # Values to represent four quadrants of a symmetric beta distribution
    categories := fnDiscretizeBeta(0.2, 0.2, 4)
    print(categories)
## references

- citation: Majumder & Bhattacharjee. 1973. Algorithm AS63. Applied Statistics, 22.
  doi: NULL
  url: NULL

- citation: WRIGHT, A. M., LLOYD, G. T. and HILLIS, D. H. 2016. Modeling character change heterogeneity in phylogenetic analyses of morphology through the use of priors. _Systematic Biology_, 65, 602–11.
  doi: 10.1093/sysbio/syv122
  url: https://doi.org/10.1093/sysbio/syv122
