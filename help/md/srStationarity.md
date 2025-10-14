## name
srStationarity
## title
Stationarity stopping rule
## description
Terminates an MCMC run when the difference between the means of individual runs
and the mean of the sample pooled from all runs ceases to be significant.
## details
This convergence criterion evaluates whether the mean of the sample pooled
from multiple runs lies outside the confidence interval of width 1 - `prob`
constructed for the mean of each individual run. Accordingly, it can only be
calculated when two or more independent runs are performed.

The number of samples to be removed as burnin before calculating the test
statistic is determined using the `burninMethod`. Different burnin lengths are
tested, increasing from 0 to 50% (for `ESS`) or 100% (for `SEM`) of the length
of the trace in increments of 10 samples. If the `ESS` option is chosen
(default), effective sample sizes (ESS) are calculated for all monitored
parameters after removing the number of samples corresponding to each candidate
burnin length. The best burnin length for a given parameter is the one that
maximizes its ESS value. If the `SEM` option is chosen, the standard error of
the mean (SEM) is calculated instead, and the best burnin length for a given
parameter is the one that minimizes its SEM value. In both cases, the final
burnin length is set to the maximum of the parameter-specific burnin lengths.

See also the tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/).
## authors
## see_also
srGelmanRubin
srGeweke
srMinESS
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
    
    # Stop when stationarity has been attained at confidence level gamma = 0.25
    stopping_rules[1] = srStationarity(prob = 0.25, file = paramFile, freq = 1000)
    
    # Create the MCMC object
    # Set nruns = 2 to ensure the stationarity statistic is applicable
    mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)
    
    # Begin the MCMC run
    mymcmc.run(rules = stopping_rules)

## references
- citation: Hill SD, Spall JC (2011). Stationarity and convergence of the Metropolis-Hastings algorithm: insights into theoretical aspects. IEEE Control Systems Magazine 39(1):56--67.
  doi: 10.1109/MCS.2018.2876959
