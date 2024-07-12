## name
srGelmanRubin
## title
Gelman–Rubin (PSRF) stopping rule

## description
Allow an MCMC run to terminate once the specified criterion has been met.
The Gelman–Rubin rule compares the variance between chains with the variance within chains; its value tends to unity (1) as runs converge.  It is widely referred to as the "potential scale reduction factor" (PSRF).

## details
## authors
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
- citation: Gelman, A; Rubin, D.B. (1992). Inference from Iterative Simulation Using Multiple Sequences. Statistical Science. 7 (4): 457–472
  doi: 10.1214/ss/1177011136
  url: null
- citation: Vats, D.; Knudson, C. Revisiting the Gelman–Rubin Diagnostic. Statist. Sci. 36 (4) 518 - 529, 2021.
  doi: 10.1214/20-STS812
  