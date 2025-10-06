## name
srGelmanRubin
## title
Gelman-Rubin (PSRF) stopping rule

## description
Allow an MCMC run to terminate once the specified criterion has been met.
The Gelman-Rubin rule compares the variance between runs with the variance within runs; its value tends to unity (1) as runs converge. It is widely referred to as the "potential scale reduction factor" (PSRF).

## details
Because the statistic is defined by comparing the variation between different runs to the variance within each run, it can only be calculated when multiple independent runs are performed, by setting the `nruns` argument to `mcmc` or `mcmcmc` to a value greater than one.

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

# Stop when the potential scale reduction factor falls below 1.01
stopping_rules[1] = srGelmanRubin(1.01, file = paramFile, freq = 1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```

## references
- citation: Gelman A, Rubin DB (1992). Inference from iterative simulation using multiple sequences. Statistical Science, 7(4):457--472.
  doi: 10.1214/ss/1177011136
- citation: Vats D, Knudson C (2021). Revisiting the Gelman--Rubin diagnostic. Statistical Science, 36(4):518--529.
  doi: 10.1214/20-STS812
  
