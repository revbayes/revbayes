## name
mnModel
## title
Monitor of Numeric Parameters
## description
Saves the values of numeric parameters sampled by an MCMC analysis to a file.
## details
Every `printgen` iterations, the `mnModel` monitor saves the values of numeric
parameters or vectors thereof (including nested vectors) to a file known as
a trace file, created at the path specified by the `filename` argument. Each
row of a trace file corresponds to one recorded iteration of an MCMC simulation
(or of an alternative sampler such as MCMCMC). The first column of a trace file
records the "Iteration" at which the samples were taken. By default, the next
three columns record the "Posterior", "Likelihood", and "Prior" of the
corresponding state of the simulation. Each of these can be suppressed by
setting the corresponding argument (`posterior`, `likelihood`, `prior`) to
`FALSE`. The subsequent columns represent individual variables. By default,
`mnModel` automatically monitors all non-constant (i.e., both stochastic and
deterministic) numeric variables in the graphical model. This behavior can be
suppressed by setting `stochasticOnly=TRUE`, and individual variables can be
omitted from the trace file using the `exclude` argument.

By default, `mnModel` prints tab-separated files, but it can also output JSON
files (by setting `format="json"`) or employ a different column separator
(specified using the `separator` argument) to produce tabular files in
alternative formats such as CSV. Tabular trace files produced by `mnModel` can
be read back into RevBayes using the `readTrace` function. The tab-separated
trace files produced under default settings can also be parsed by a number
of third-party tools such as Tracer (Rambaut et al. 2018) or RWTY (Warren et
al. 2017).

The `mnModel` can only save the values of simple, numeric parameters. More
complex parameters such as trees or rate matrices can instead be recorded using
the `mnFile` monitor.
 
## authors
## see_also
mnFile
readTrace
Trace
## example
    # Binomial example: estimate success probability given 7 successes out of 20 trials
    r ~ dnExp(10)
    p := Probability(ifelse(r < 1, r, 1))
    n <- 20
    k ~ dnBinomial(n, p)
    k.clamp(7)
    mymodel = model(k)

    moves = VectorMoves()
    moves.append( mvScale(r, weight=1) )

    # Set up a monitor for both r and p
    all_params = mnModel(filename="all_params.log", printgen=10)
    # Set up a monitor for r only
    stoch_only = mnModel(filename="stoch_only.log", stochasticOnly=TRUE, printgen=10)
    # Set up a monitor for p only
    p_only = mnModel(filename="p_only.log", exclude=["r"], printgen=10)

    # Apply the monitors and run a short simulation
    mymcmc = mcmc(model=mymodel, moves=moves, monitors=[all_params, stoch_only, p_only])
    mymcmc.run(generations=1000)

## references
- citation: Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA (2018). Posterior summarization in Bayesian phylogenetics using Tracer 1.7. Systematic Biology, 67(5):901-904.
  doi: 10.1093/sysbio/syy032
  url: https://academic.oup.com/sysbio/article/67/5/901/4989127
- citation: Warren DL, Geneva AJ, Lanfear R (2017). RWTY (R We There Yet): an R package for examining convergence of Bayesian phylogenetic analyses. Molecular Biology and Evolution, 34(4):1016-1020.
  doi: 10.1093/molbev/msw279
  url: https://academic.oup.com/mbe/article/34/4/1016/2900564
