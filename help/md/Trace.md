## name
Trace
## title
Trace of numeric parameter values.
## description
Corresponds to a single column of a log file, usually produced by the `mnModel` monitor in an MCMC or MCMCMC run.
## details
Method description:

- `getBurnin()`: Return the number of samples that were discarded as burnin, after thinning.
- `getValues()`: Return a vector (e.g. `Real[]`) containing the values of the sampled parameter,
  after excluding burnin samples.
- `setBurnin()`: Modify the number (if >= 1) or fraction (if < 1) of samples to discard as burnin, after thinning.
- `size()`, `getNumberSamples()`: Report the number of values stored in `trace`,
  including (`post = FALSE`) or excluding (`post = TRUE`) burnin samples.
- `summarize()`: Display summary statistics of trace.

## authors
## see_also
mnModel
readTrace
TraceTree
## example
    # Read in a log file as a vector of traces
    traces = readTrace("out.log", burnin = 0)
    
    # Get the posterior trace (2nd column in the log file)
    posterior = traces[2]
    posterior.getBurnin() # will return 0
    
    # Change burnin from 0 to 50 samples
    posterior.setBurnin(50)
    posterior.getBurnin() # will return 50
    
    # Get summary statistics
    posterior.summarize()

## references
