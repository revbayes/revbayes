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
# Specify a model
mymodel = model(tree)

# Define monitors to monitor parameters for convergence
paramFile = "parameters.log"
# Configure the parameter file to monitor only the parameters whose convergence shall be assessed
monitors[1] = mnModel( filename = paramFile, printgen = 10, stochasticOnly = TRUE, exclude = ["rel_br_lengths"] )

# Define threshold rules (optional): Stop as soon as ANY has been met
stopping_rules[1] = srMaxIteration(200000)
stopping_rules[2] = srMaxTime(15, "hours")

# Check for convergence every 1000 iterations
convergenceFreq = 1000

# Add convergence rules: Stop as soon as ALL have been met
stopping_rules[3] = srMinESS(50, file = paramFile, freq = convergenceFreq)
stopping_rules[4] = srGelmanRubin(1.01, file = paramFile, freq = convergenceFreq)
stopping_rules[5] = srGeweke(prob = 0.001, file = paramFile, freq = convergenceFreq)
stopping_rules[6] = srStationarity(prob = 0.01, file= paramFile ,freq = convergenceFreq)

mymcmc = mcmc(mymodel, monitors, moves)
# Begin the MCMC run.
# The checkpoint file will allow the run to be resumed from where it left off
# if you later opt for more stringent stopping rules
mymcmc.run(rules = stopping_rules, checkpointFile = "checkpoint.ckp", checkpointInterval = 1000)
```

## references
- citation: Geweke, J. Evaluating the accuracy of sampling-based approaches to calculating posterior moments.
\ In Bayesian Statistics 4 (ed JM Bernado, JO Berger, AP Dawid and AFM Smith). 
\ Clarendon Press, Oxford, UK.
