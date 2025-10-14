## name
srGelmanRubin
## title
Gelman–Rubin (PSRF) stopping rule
## description
Terminates an MCMC run when the Gelman-Rubin statistic drops below the
specified value.
## details
The Gelman–Rubin statistic, also referred to as the potential scale reduction
factor (PSRF), compares the variance of the sample pooled from multiple runs to
the sum of variances calculated from individual runs. Accordingly, it can only
be calculated when two or more independent runs are performed, and its value
tends to unity (1) as the runs converge.

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
srGeweke
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
    
    # Stop when the potential scale reduction factor falls below 1.01
    stopping_rules[1] = srGelmanRubin(1.01, file = paramFile, freq = 1000)
    
    # Create the MCMC object.
    # Set nruns = 2 to ensure the Gelman-Rubin statistic is applicable
    mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)
    
    # Begin the MCMC run
    mymcmc.run(rules = stopping_rules)

## references
- citation: Gelman A, Rubin DB (1992). Inference from iterative simulation using multiple sequences. Statistical Science, 7(4):457--472.
  doi: 10.1214/ss/1177011136
- citation: Vats D, Knudson C (2021). Revisiting the Gelman--Rubin diagnostic. Statistical Science, 36(4):518--529.
  doi: 10.1214/20-STS812
  
