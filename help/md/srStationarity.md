## name
srStationarity
## title

## description
Allow an MCMC run to terminate once the specified criterion has been met.

## details
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

# Stop when stationarity has been accomplished at significance level = 0.05
stopping_rules[1] = srStationarity(prob = 0.05, file = paramFile, freq = 1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```

## references
