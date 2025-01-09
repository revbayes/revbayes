## name
srStationarity
## title
Stationarity stopping rule

## description
Allow an MCMC run to terminate once the specified criterion has been met.
An MCMC sample can be considered stationary once its mean, variance and autocorrelation structure do not change over time.

## details
Because the statistic is defined by comparing different runs, it can only be calculated when multiple independent runs are performed, by setting the `nruns` argument to `mcmc` or `mcmcmc` to a value greater than one.

## authors
## see_also
- Tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/)

## example
```
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
mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```

## references
- citation: Hill, S.D. and Spall, J.C. 2011. Stationarity and Convergence of the Metropolis-Hastings Algorithm: Insights into Theoretical Aspects. IEEE Control Systems Magazine 39.
  doi: 10.1109/MCS.2018.2876959
  url: null
