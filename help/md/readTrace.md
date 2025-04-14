## name
readTrace
## title
Read an MCMC log file.
## description
Reads parameter values from a log file,  usually produced as the output of an MCMC or MCMCMC run.
## details

Read an MCMC log file with field delimited by `separator`.
Then drop the first `burnin` iterations if `burnin` is an integer,
or the fraction `burnin` of iterations if `burnin` if a Real number.
Then we keep every $n$th entry if the `thinning` is $n$.

Pertinent methods of a `Trace` object include:

- `trace.getBurnin()`: Return the number of samples that were discarded as burnin, after thinning.
- `trace.getValues()`: Return a vector (e.g. `Real[]`) containing the values of the sampled parameter,
  after excluding burnin samples.
- `trace.setBurnin()`: Intended to allow modification of the burnin size. Currently has no effect.
- `trace.size()`, `trace.getNumberSamples()`: Report the number of values stored in `trace`,
  including (`post = FALSE`) or excluding (`post = TRUE`) burnin samples.
- `trace.summarize()`: Display summary statistics of trace.

## authors
## see_also
readTreeTrace
readCharacterData
## example
trace <- readTrace(filename, burnin=burnin)[1]
## references
