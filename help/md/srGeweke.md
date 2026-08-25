## name
srGeweke
## title
Geweke stopping rule
## description
Terminates an MCMC run when the Geweke test statistic ceases to be significant
at the specified level.
## details
Geweke (1992) proposed a convergence diagnostic for Markov chains based on
a test for equality of the means of the first and last part of a Markov chain
(by default the first 10% and the last 50%). If the samples are drawn from the
stationary distribution of the chain, the two means are equal and Geweke's
statistic has an asymptotically standard normal distribution.

The test statistic is a standard Z-score: the difference between the two sample
means divided by its estimated standard error. The standard error is estimated
from the spectral density at zero and so accounts for any autocorrelation. The
Z-score is calculated under the assumption that the two parts of the chain are
asymptotically independent, which requires that the sum of `frac1` and `frac2`
be strictly less than 1.

The number of samples to be removed as burnin before calculating the test
statistic is determined using the `burninMethod`. If the `"ESS"` (default) or
`"SEM"` options are chosen, different burnin lengths are tested, increasing
from 0 to 50% (for `"ESS"`) or 100% (for `"SEM"`) of the length of the trace in
increments of 10 samples. The `"ESS"` option calculates effective sample sizes
(ESS) for all monitored parameters after removing the number of samples
corresponding to each candidate burnin length. The best burnin length for a
given parameter is the one that maximizes its ESS value. The `"SEM"` option
instead calculates the standard error of the mean (SEM), and the best burnin
length for a given parameter is the one that minimizes its SEM value. In both
cases, the final burnin length is set to the maximum of the parameter-specific
burnin lengths.

Alternatively, the user may set `burninMethod` to `"fixed"`, which discards a
constant fraction of the samples collected up to that point. This fraction can
be specified using the `burnin` argument, and is set to 0.25 by default. The
`burnin` argument has no effect if the `"ESS"` or `"SEM"` options are chosen,
and a corresponding warning is displayed if the user explicitly sets the
argument without specifying `burninMethod="fixed"`. The `"fixed"` option is 
appropriate for analyses with very long parameter traces and large numbers of
monitored variables, for which the automatic burnin determination may be too
computationally demanding.

See also the tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/).

This help file incorporates text by Martyn Plummer.
## authors
## see_also
mcmc
srGelmanRubin
srMinESS
srStationarity
## example
    # Binomial example: estimate success probability given 7 successes out of 20 trials
    r ~ dnExp(10)
    p := Probability(ifelse(r < 1, r, 1))
    n <- 20
    k ~ dnBinomial(n, p)
    k.clamp(7)
    mymodel = model(k)
    
    moves = VectorMoves()
    moves.append( mvSlide(r, delta=0.1, weight=1) )
    
    paramFile = "parameters.log"
    
    monitors = VectorMonitors()
    monitors.append( mnModel(filename=paramFile, printgen=100, p) )
    
    # Stop when the Geweke test statistic is no longer significant at alpha = 0.001
    stopping_rules[1] = srGeweke( prob=0.001, file=paramFile, freq=10000 )
    
    # Create the MCMC object
    mymcmc = mcmc( mymodel, monitors, moves )
    
    # Begin the MCMC run
    mymcmc.run( rules = stopping_rules )

## references
- citation: Geweke J (1992). Evaluating the accuracy of sampling-based approaches to calculating posterior moments. Pp. 169--194 in Bernado M, Berger JO, Dawid AP, Smith AFM, eds. Bayesian Statistics 4. Clarendon Press, Oxford, UK.
  doi: 10.1093/oso/9780198522669.003.0010
  url: https://academic.oup.com/book/54041/chapter-abstract/422209572
