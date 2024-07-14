## name
srGeweke
## title
Geweke stopping rule

## description
Allow an MCMC run to terminate once the specified criterion has been met.

Geweke (1992) proposed a convergence diagnostic for Markov chains based on a test for equality of the means of the first and last part of a Markov chain (by default the first 10% and the last 50%). If the samples are drawn from the stationary distribution of the chain, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution.

The test statistic is a standard Z-score: the difference between the two sample means divided by its estimated standard error. The standard error is estimated from the spectral density at zero and so takes into account any autocorrelation.

The Z-score is calculated under the assumption that the two parts of the chain are asymptotically independent, which requires that the sum of `frac1` and `frac2` be strictly less than 1.

## details
## authors
Incorporates text by Martyn Plummer

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

# Stop when the Geweke test statistic becomes significant at alpha = 0.001
stopping_rules[1] = srGeweke( prob=0.001, file=paramFile, freq=10000 )

# Create the MCMC object
mymcmc = mcmc( mymodel, monitors, moves )

# Begin the MCMC run
mymcmc.run( rules = stopping_rules )
```

## references
- citation: Geweke, J. Evaluating the accuracy of sampling-based approaches to calculating posterior moments.
\ In Bayesian Statistics 4 (ed JM Bernado, JO Berger, AP Dawid and AFM Smith). 
\ Clarendon Press, Oxford, UK.
