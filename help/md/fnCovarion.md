## name
fnCovarion
## title
The Covarion model rate matrix.
## description
The `fnCovarion` function defines a covarion model rate matrix for character evolution.
The resulting rate matrix incorporates rate heterogeneity where the characters/sites are allowed to move between rate categories with a switching rate.

## details
The covarion model allows for variation in evolutionary rates across sites over time, accommodating shifts in character state evolution.

## authors
Sebastian Hoehna and Lyndon Coghill

## see_also
fnJC

fnF81

## example
    # define number of rate categories
    num_cats = 2

    # define rate scalars
    rate_scalars <- [ 0, 1.0 ]

    # obtain substitution model rate matrix for rescaling
    for ( i in 1:num_cats ) {
      Q_sub[i] := fnJC( 4 )
    }

    # define switching rate
    switching_rate = 0.1

    for ( i in 1:num_cats ) {
      for ( j in 1:num_cats ) {
        if ( i == j ) {
          switch_rates[i][j] := 0
        } else {
          switch_rates[i][j] := switching_rate
        }
      }
    }

    # finally construct the covarion rate matrix
    Q_Cov := fnCovarion(RateMatrices = Q_sub, RateScalars = rate_scalars, SwitchRates = switch_rates, rescaled = FALSE)


## references
- citation: Fitch, W. M., & Markowitz, E. (1970). An improved method for determining codon variability in a gene and its application to the rate of fixation of mutations in evolution. _Biochemical genetics_, 4, 579-593.
  doi: 10.1007/BF00486096
  url: https://doi.org/10.1007/BF00486096

- citation: Tuffley, C., & Steel, M. (1998). Modeling the covarion hypothesis of nucleotide substitution. _Mathematical biosciences_, 147(1), 63-91.
  doi: 10.1016/S0025-5564(97)00081-3
  url: https://doi.org/10.1016/S0025-5564(97)00081-3
